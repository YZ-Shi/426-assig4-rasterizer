"use strict";
var Reflection = Reflection || {
  ambient: new Pixel(0, 0, 0),
  diffuse: new Pixel(1.0, 1.0, 1.0),
  specular: new Pixel(1.0, 1.0, 1.0),
  shininess: 20,
};

Reflection.phongReflectionModel = function(vertex, view, normal, lightPos, phongMaterial) {
  var color = new Pixel(0, 0, 0);
  normal.normalize();

  // diffuse
  var light_dir = new THREE.Vector3().subVectors(lightPos, vertex).normalize();
  var ndotl = normal.dot(light_dir);
  color.plus(phongMaterial.diffuse.copy().multipliedBy(ndotl));

  // Ambient color and specular color

  color.plus(phongMaterial.ambient); // ambient color
  //specular color
  let r = light_dir.clone().negate().reflect(normal);
  let v = view.clone().sub(vertex).normalize();
  let rdotv = Math.max(0, r.dot(v));
  color.plus(phongMaterial.specular.copy().multipliedBy(rdotv ** phongMaterial.shininess));


  return color;
};

var Renderer = Renderer || {
  meshInstances: new Set(),
  width: 1000,
  height: 750,
  negNear: 0.3,
  negFar: 1000,
  fov: 45,
  lightPos: new THREE.Vector3(10, 10, -10),
  shaderMode: "",
  cameraLookAtVector: new THREE.Vector3(0, 0, 0),
  cameraPosition: new THREE.Vector3(0, 0, -10),
  cameraUpVector: new THREE.Vector3(0, -1, 0),
  cameraUpdated: true,
};

Renderer.updateCameraParameters = function() {
  this.camera.position.copy(this.cameraPosition);
  this.camera.up.copy(this.cameraUpVector);
  this.camera.lookAt(this.cameraLookAtVector);
};

Renderer.initialize = function() {
  this.buffer = new Image(this.width, this.height);
  this.zBuffer = [];

  // set camera
  this.camera = new THREE.PerspectiveCamera(
    this.fov,
    this.width / this.height,
    this.negNear,
    this.negFar
  );
  this.updateCameraParameters();

  this.clearZBuffer();
  this.buffer.display(); // initialize canvas
};

Renderer.clearZBuffer = function() {
  for (var x = 0; x < this.width; x++) {
    this.zBuffer[x] = new Float32Array(this.height);
    for (var y = 0; y < this.height; y++) {
      this.zBuffer[x][y] = 1; // z value is in [-1 1];
    }
  }
};

Renderer.addMeshInstance = function(meshInstance) {
  assert(meshInstance.mesh, "meshInstance must have mesh to be added to renderer");
  this.meshInstances.add(meshInstance);
};

Renderer.removeMeshInstance = function(meshInstance) {
  this.meshInstances.delete(meshInstance);
};

Renderer.clear = function() {
  this.buffer.clear();
  this.clearZBuffer();
  Main.context.clearRect(0, 0, Main.canvas.width, Main.canvas.height);
};

Renderer.displayImage = function() {
  this.buffer.display();
};

Renderer.render = function() {
  this.clear();

  var eps = 0.01;
  if (
    !(
      this.cameraUpVector.distanceTo(this.camera.up) < eps &&
      this.cameraPosition.distanceTo(this.camera.position) < eps &&
      this.cameraLookAtVector.distanceTo(Main.controls.target) < eps
    )
  ) {
    this.cameraUpdated = false;
    // update camera position
    this.cameraLookAtVector.copy(Main.controls.target);
    this.cameraPosition.copy(this.camera.position);
    this.cameraUpVector.copy(this.camera.up);
  } else {
    // camera's stable, update url once
    if (!this.cameraUpdated) {
      Gui.updateUrl();
      this.cameraUpdated = true; //update one time
    }
  }

  this.camera.updateMatrixWorld();
  this.camera.matrixWorldInverse.getInverse(this.camera.matrixWorld);

  // light goes with the camera, COMMENT this line for debugging if you want
  this.lightPos = this.camera.position;

  for (var meshInst of this.meshInstances) {
    var mesh = meshInst.mesh;
    if (mesh !== undefined) {
      for (var faceIdx = 0; faceIdx < mesh.faces.length; faceIdx++) {
        var face = mesh.faces[faceIdx];
        var verts = [mesh.vertices[face.a], mesh.vertices[face.b], mesh.vertices[face.c]];
        var vert_normals = [
          mesh.vertex_normals[face.a],
          mesh.vertex_normals[face.b],
          mesh.vertex_normals[face.c],
        ];

        // camera's view matrix = K * [R | t] where K is the projection matrix and [R | t] is the inverse of the camera pose
        var viewMat = new THREE.Matrix4().multiplyMatrices(
          this.camera.projectionMatrix,
          this.camera.matrixWorldInverse
        );

        Renderer.drawTriangle(verts, vert_normals, mesh.uvs[faceIdx], meshInst.material, viewMat);
      }
    }
  }

  this.displayImage();
};

Renderer.getPhongMaterial = function(uv_here, material) {
  var phongMaterial = {};
  phongMaterial.ambient = Reflection.ambient;

  if (material.diffuse === undefined || uv_here === undefined) {
    phongMaterial.diffuse = Reflection.diffuse;
  } else if (Pixel.prototype.isPrototypeOf(material.diffuse)) {
    phongMaterial.diffuse = material.diffuse;
  } else {
    // note that this function uses point sampling. it would be better to use bilinear
    // subsampling and mipmaps for area sampling, but this good enough for now...
    phongMaterial.diffuse = material.diffuse.getPixel(
      Math.floor(uv_here.x * material.diffuse.width),
      Math.floor(uv_here.y * material.diffuse.height)
    );
  }

  if (material.specular === undefined || uv_here === undefined) {
    phongMaterial.specular = Reflection.specular;
  } else if (Pixel.prototype.isPrototypeOf(material.specular)) {
    phongMaterial.specular = material.specular;
  } else {
    phongMaterial.specular = material.specular.getPixel(
      Math.floor(uv_here.x * material.specular.width),
      Math.floor(uv_here.y * material.specular.height)
    );
  }

  phongMaterial.shininess = Reflection.shininess;

  return phongMaterial;
};

Renderer.projectVerticesNaive = function(verts) {
  // this is a naive orthogonal projection that does not even consider camera pose
  var projectedVerts = [];

  var orthogonalScale = 5;
  for (var i = 0; i < 3; i++) {
    projectedVerts[i] = new THREE.Vector4(verts[i].x, verts[i].y, verts[i].z, 1.0);

    projectedVerts[i].x /= orthogonalScale;
    projectedVerts[i].y /= (orthogonalScale * this.height) / this.width;

    projectedVerts[i].x = (projectedVerts[i].x * this.width) / 2 + this.width / 2;
    projectedVerts[i].y = (projectedVerts[i].y * this.height) / 2 + this.height / 2;
  }

  return projectedVerts;
};

Renderer.projectVertices = function(verts, viewMat) {
  // Vector3/Vector4 array of projected vertices in screen space coordinates
  // (you still need z for z buffering)
  var projectedVerts = [];


  for (var i = 0; i < 3; i++) {
    projectedVerts[i] = new THREE.Vector4(verts[i].x, verts[i].y, verts[i].z, 1.0);

    projectedVerts[i].applyMatrix4(viewMat);
    projectedVerts[i].divideScalar(projectedVerts[i].w);

    if (projectedVerts[i].z < this.negNear || projectedVerts[i].z > this.negFar)
    return undefined;

    projectedVerts[i].x = (projectedVerts[i].x + 1) / 2 * this.width;
    projectedVerts[i].y = (projectedVerts[i].y + 1) / 2 * this.height;
  }


  return projectedVerts;
};

Renderer.computeBoundingBox = function(projectedVerts) {
  // Compute the screen-space bounding box for the triangle defined in projectedVerts[0-2].
  // We will need to call this helper function in the shading functions
  // to loop over pixel locations in the bounding box for rasterization.

  var box = {};
  box.minX = -1;
  box.minY = -1;
  box.maxX = -1;
  box.maxY = -1;


  let minX = Math.min(projectedVerts[0].x, projectedVerts[1].x, projectedVerts[2].x);
  let minY = Math.min(projectedVerts[0].y, projectedVerts[1].y, projectedVerts[2].y);
  let maxX = Math.max(projectedVerts[0].x, projectedVerts[1].x, projectedVerts[2].x);
  let maxY = Math.max(projectedVerts[0].y, projectedVerts[1].y, projectedVerts[2].y);

  box.minX = Math.max(0, Math.floor(minX));
  box.minY = Math.max(0, Math.floor(minY));
  box.maxX = Math.min(this.width - 1, Math.ceil(maxX));
  box.maxY = Math.min(this.height - 1, Math.ceil(maxY));

  return box;
};

Renderer.computeBarycentric = function(projectedVerts, x, y) {
  var triCoords = [];
  // (see https://fgiesen.wordpress.com/2013/02/06/the-barycentric-conspirac/)
  // return undefined if (x,y) is outside the triangle

  let wa = ((projectedVerts[1].x - projectedVerts[2].x) * (projectedVerts[2].y - y)
  - (projectedVerts[2].x - x) * (projectedVerts[1].y - projectedVerts[2].y))
  / ((projectedVerts[0].x - projectedVerts[2].x) * (projectedVerts[1].y - projectedVerts[2].y)
  - (projectedVerts[1].x - projectedVerts[2].x) * (projectedVerts[0].y - projectedVerts[2].y));
  if (wa < 0) return undefined;
  let wb = ((projectedVerts[0].x - projectedVerts[2].x) * (projectedVerts[2].y - y)
  - (projectedVerts[2].x - x) * (projectedVerts[0].y - projectedVerts[2].y))
  / ((projectedVerts[1].x - projectedVerts[2].x) * (projectedVerts[0].y - projectedVerts[2].y)
  - (projectedVerts[0].x - projectedVerts[2].x) * (projectedVerts[1].y - projectedVerts[2].y));
  if (wb < 0) return undefined;
  let wc = 1 - wa - wb;
  if (wc < 0) return undefined;
  triCoords = [wa, wb, wc];

  return triCoords;
};

Renderer.drawTriangleWire = function(projectedVerts) {
  var color = new Pixel(1.0, 0, 0);
  for (var i = 0; i < 3; i++) {
    var va = projectedVerts[(i + 1) % 3];
    var vb = projectedVerts[(i + 2) % 3];

    var ba = new THREE.Vector2(vb.x - va.x, vb.y - va.y);
    var len_ab = ba.length();
    ba.normalize();
    // draw line
    for (var j = 0; j < len_ab; j += 0.5) {
      var x = Math.round(va.x + ba.x * j);
      var y = Math.round(va.y + ba.y * j);
      this.buffer.setPixel(x, y, color);
    }
  }
};

Renderer.drawTriangleFlat = function(verts, projectedVerts, normals, uvs, material) {
  // Flat shader
  // Color of each face is computed based on the face normal
  // (average of vertex normals) and face centroid.

  // compute face normal and centroid
  let normal = new THREE.Vector3().addVectors(normals[0], normals[1]).add(normals[2]).divideScalar(3).normalize();
  let centroid = new THREE.Vector3().addVectors(verts[0], verts[1]).add(verts[2]).divideScalar(3);
  let phongMaterial;
  // texture mapping if applicable
  if (uvs !== undefined) {
    let uv = {x: 0, y: 0};
    let centroidProj = new THREE.Vector3().addVectors(projectedVerts[0], projectedVerts[1]).add(projectedVerts[2]).divideScalar(3);
    let bCoords = this.computeBarycentric(projectedVerts, centroidProj.x, centroidProj.y);
    uv.x = bCoords[0] * uvs[0].x + bCoords[1] * uvs[1].x + bCoords[2] * uvs[2].x;
    uv.y = bCoords[0] * uvs[0].y + bCoords[1] * uvs[1].y + bCoords[2] * uvs[2].y;
    phongMaterial = this.getPhongMaterial(uv, material);
  }
  else phongMaterial = this.getPhongMaterial(uvs, material);
  // derive color from phong reflection model
  let color = Reflection.phongReflectionModel(centroid, this.cameraPosition, normal, this.lightPos, phongMaterial);
  // iterate through bounding box and render pixels if within the triangel and closer
  let box = this.computeBoundingBox(projectedVerts);
  for (let x = box.minX; x <= box.maxX; x++) {
    for (let y = box.minY; y <= box.maxY; y++) {
      let bCoords = this.computeBarycentric(projectedVerts, x, y);
      if (bCoords === undefined) continue; // out of the triangle
      // Check if z is closer
      let z = bCoords[0] * projectedVerts[0].z + bCoords[1] * projectedVerts[1].z + bCoords[2] * projectedVerts[2].z;
      if (z >= this.zBuffer[x][y]) continue;
      this.buffer.setPixel(x, y, color);
      this.zBuffer[x][y] = z;
    }
  }

};

Renderer.drawTriangleGouraud = function(verts, projectedVerts, normals, uvs, material) {
  // Gouraud shader
  // Interpolate the color for each pixel in the triangle using the barycentric coordinate.

  // find colors at the vertices
  let colors = [];
  if (uvs === undefined) {
    let phongMaterial = this.getPhongMaterial(uvs, material);
    for (let i = 0; i < 3; i++) {
      colors[i] = Reflection.phongReflectionModel(verts[i], this.cameraPosition, normals[i], this.lightPos, phongMaterial);
    }
  } else { // texture mapping
    for (let i = 0; i < 3; i++) {
      let phongMaterial = this.getPhongMaterial(uvs[i], material);
      colors[i] = Reflection.phongReflectionModel(verts[i], this.cameraPosition, normals[i], this.lightPos, phongMaterial);
    }
  }
  let box = this.computeBoundingBox(projectedVerts);
  for (let x = box.minX; x <= box.maxX; x++) {
    for (let y = box.minY; y <= box.maxY; y++) {
      let bCoords = this.computeBarycentric(projectedVerts, x, y);
      if (bCoords === undefined) continue; // out of the triangle
      // Check if z is closer
      let z = bCoords[0] * projectedVerts[0].z + bCoords[1] * projectedVerts[1].z + bCoords[2] * projectedVerts[2].z;
      if (z >= this.zBuffer[x][y]) continue;
      // find color
      let color = new Pixel(0, 0, 0);
      for (let j = 0; j < 3; j++) {
        color.plus(colors[j].copyMultiplyScalar(bCoords[j]));
      }
      this.buffer.setPixel(x, y, color);
      this.zBuffer[x][y] = z;
    }
  }

};

Renderer.drawTrianglePhong = function(verts, projectedVerts, normals, uvs, material) {
  // Phong shader
  // (1) Basic Phong shader: Interpolate the normal and vertex for each pixel in the triangle
  //                         using the barycentric coordinate.
  // (2) Texture mapping: If uvs is provided, compute interpolated uv coordinates
  //                      and map the phong material texture (if available)
  //                      at the uv coordinates to the pixel location.
  // (3) XYZ normal mapping: If xyz normal texture exists for the material,
  //                         convert the RGB value of the XYZ normal texture at the uv coordinates
  //                         to a normal vector and apply it at the pixel location.

  let phongMaterial;
  if (uvs === undefined) { // texture mapping
    phongMaterial = this.getPhongMaterial(uvs, material);
  }
  let box = this.computeBoundingBox(projectedVerts);
  let flag;
  for (let y = box.minY; y <= box.maxY; y++) {
    for (let x = box.minX; x <= box.maxX; x++) {
      let bCoords = this.computeBarycentric(projectedVerts, x, y);
      if (bCoords === undefined) { // out of the triangle
        if (flag) { // reached the right end of the triangle
          flag = false;
          break; // reaches right side of triangle
        }
        continue;
      }
      flag = true;
      // Check if z is closer
      let z = bCoords[0] * projectedVerts[0].z + bCoords[1] * projectedVerts[1].z + bCoords[2] * projectedVerts[2].z;
      if (z >= this.zBuffer[x][y]) continue;
      // find normal and positions
      let normal = new THREE.Vector3();
      let position = new THREE.Vector3();
      let uv = {x: 0, y: 0};
      for (let j = 0; j < 3; j++) {
        normal.add(normals[j].clone().multiplyScalar(bCoords[j]));
        position.add(verts[j].clone().multiplyScalar(bCoords[j]));
        if (uvs !== undefined) {
          uv.x += uvs[j].x * bCoords[j];
          uv.y += uvs[j].y * bCoords[j];
        }
      }
      // normal mapping
      let normXyz;
      if (material.xyzNormal !== undefined) {
        let w = Math.floor(uv.x * material.xyzNormal.width);
        let h = Math.floor(uv.y * material.xyzNormal.height);
        let rgb = material.xyzNormal.getPixel(w, h);
        let rgbMapped = new Pixel(-1, -1, -1);
        rgbMapped.plus(rgb.copyMultiplyScalar(2));
        normXyz = new THREE.Vector3(rgbMapped.r, rgbMapped.g, rgbMapped.b);
        normXyz.normalize();
      }
      if (uvs !== undefined) {
        phongMaterial = this.getPhongMaterial(uv, material);
      }
      let color;
      if (material.xyzNormal !== undefined) {
        color = Reflection.phongReflectionModel(position, this.cameraPosition, normXyz, this.lightPos, phongMaterial);
      }
      else color = Reflection.phongReflectionModel(position, this.cameraPosition, normal, this.lightPos, phongMaterial);
      // find color
      this.buffer.setPixel(x, y, color);
      this.zBuffer[x][y] = z;
    }
  }
  // sort the vertices based on the y coordinate
  /*let projectedSort = [];
  for (let i = 0; i < 3; i++) {
    projectedSort[i] = projectedVerts[i];
  }
  projectedSort.sort(function(a, b){return a.y-b.y});
  // calculate the slopes of each of the edges
  let delta01 = (projectedSort[1].x - projectedSort[0].x) / (projectedSort[1].y - projectedSort[0].y);
  let delta02 = (projectedSort[2].x - projectedSort[0].x) / (projectedSort[2].y - projectedSort[0].y);
  let delta12 = (projectedSort[2].x - projectedSort[1].x) / (projectedSort[2].y - projectedSort[1].y);
  let xL, xR, deltaL, deltaR;
  xL = Math.floor(projectedSort[0].x);
  xR = Math.ceil(projectedSort[0].x);
  if (projectedSort[0].y - projectedSort[1].y < 1) { // flat top edge
    if (projectedSort[0].x < projectedSort[1].x) xR = projectedSort[1].x;
    else xL = projectedSort[1].x;
  }

  // for each scan line
  for (let y = Math.floor(projectedSort[0].y); y <= Math.ceil(projectedSort[2].y); y++) {
    if (projectedSort[0].y - projectedSort[1].y < 1) { // flat top edge
      if (projectedSort[0].x < projectedSort[1].x) {
        deltaL = delta02;
        deltaR = delta12;
      } else {
        deltaL = delta12;
        deltaR = delta02;
      }
    }
    else if (y < projectedSort[1].y) {
      if (projectedSort[1].x < projectedSort[2].x) {
        deltaL = delta01;
        deltaR = delta02;
      }
      else {
        deltaL = delta02;
        deltaR = delta01;
      }
    } else {
      if (projectedSort[1].x < projectedSort[2].x) deltaL = delta12;
      else deltaR = delta12;

    }
    // from left to right in scan line
    let flag;
    for (let x = xL; x <= box.maxX; x++) {
      let bCoords = this.computeBarycentric(projectedVerts, x, y);
      if (bCoords === undefined) {
        if (flag) break;
        continue;
      }
      flag = true;
      // Check if z is closer
      let z = bCoords[0] * projectedVerts[0].z + bCoords[1] * projectedVerts[1].z + bCoords[2] * projectedVerts[2].z;
      if (z >= this.zBuffer[x][y]) continue;
      // find normal and positions
      let normal = new THREE.Vector3();
      let position = new THREE.Vector3();
      let uv = {x: 0, y: 0};
      for (let j = 0; j < 3; j++) {
        normal.add(normals[j].clone().multiplyScalar(bCoords[j]));
        position.add(verts[j].clone().multiplyScalar(bCoords[j]));
        if (uvs !== undefined) {
          uv.x += uvs[j].x * bCoords[j];
          uv.y += uvs[j].y * bCoords[j];
        }
      }
      // normal mapping
      let normXyz;
      if (material.xyzNormal !== undefined) {
        let w = Math.floor(uv.x * material.xyzNormal.width);
        let h = Math.floor(uv.y * material.xyzNormal.height);
        let rgb = material.xyzNormal.getPixel(w, h);
        let rgbMapped = new Pixel(-1, -1, -1);
        rgbMapped.plus(rgb.copyMultiplyScalar(2));
        normXyz = new THREE.Vector3(rgbMapped.r, rgbMapped.g, rgbMapped.b);
        normXyz.normalize();
      }
      if (uvs !== undefined) {
        phongMaterial = this.getPhongMaterial(uv, material);
      }
      let color;
      if (material.xyzNormal !== undefined) {
        color = Reflection.phongReflectionModel(position, this.cameraPosition, normXyz, this.lightPos, phongMaterial);
      }
      else color = Reflection.phongReflectionModel(position, this.cameraPosition, normal, this.lightPos, phongMaterial);
      // find color
      this.buffer.setPixel(x, y, color);
      this.zBuffer[x][y] = z;
    }
    // update the left and right ends for the next sweepline
    xL = Math.floor(xL + deltaL);
    xR = Math.ceil(xR + deltaR);
  }*/

};

Renderer.drawTriangle = function(verts, normals, uvs, material, viewMat) {
  var projectedVerts = this.projectVertices(verts, viewMat);
  if (projectedVerts === undefined) {
    // not within near and far plane
    return;
  } else if (projectedVerts.length <= 0) {
    projectedVerts = this.projectVerticesNaive(verts);
  }

  switch (this.shaderMode) {
    case "Wire":
    this.drawTriangleWire(projectedVerts);
    break;
    case "Flat":
    this.drawTriangleFlat(verts, projectedVerts, normals, uvs, material);
    break;
    case "Gouraud":
    this.drawTriangleGouraud(verts, projectedVerts, normals, uvs, material);
    break;
    case "Phong":
    this.drawTrianglePhong(verts, projectedVerts, normals, uvs, material);
    break;
    default:
  }
};
