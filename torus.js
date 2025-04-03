class Point {
    point;
    dim;
    constructor(point) {
        if (point == null) {
            return;
        }
        this.point = point;
        this.dim = point.length;
    }

    getPoint() {
        return this.point;
    }

    sub(that) {
        return new Vector(new Point(that.point), new Point(this.point), null);
    }

    getDim() {
        return this.dim;
    }

    dist(that) {
        let sum = 0;
        for (let i = 0; i < this.dim; i++) {
            sum += (this.point[i] - that.point[i]) ** 2;
        }
        return Math.sqrt(sum);
    }
}

class Vector {
    start;
    end;
    coords;
    dim;
    // let ray = new Vector(camPnt.getEnd(), new Point([x, y, 0]), null);

    constructor(start, end, coords) {
        if (start != null && end != null) {
            this.dim = start.getDim();
            this.start = start;
            this.end = end;
            coords = [];
            for (let i = 0; i < this.dim; i++) {
                coords[i] = end.getPoint()[i] - start.getPoint()[i];
            }
            this.coords = coords;
            return;
        }
        if (end != null) {
            start = [];
            coords = [];
            for (let i = 0; i < end.getDim(); i++) {
                start[i] = 0;
                coords[i] = end.getPoint()[i];
            }
            this.start = new Point(start);
            this.end = end;
            this.coords = coords;
            this.dim = coords.length;
            return;
            
        }
        start = [];
        for (let i = 0; i < coords.length; i++) {
            start[i] = 0;
        }
        this.start = new Point(start);
        this.end = new Point(coords);
        this.coords = coords;
        this.dim = coords.length;
    }

    getCoords() {
        return this.coords;
    }

    getStart() {
        return this.start;
    }

    getEnd() {
        return this.end;
    }

    add(that) {
        let newEnd = [];
        for (let i = 0; i < this.dim; i++) {
            newEnd[i] = this.end.getPoint()[i] + that.coords[i];
        }
        return new Vector(this.start, new Point(newEnd), null);
    }

    sub(that) {
        let newEnd = [];
        that = that.scale(-1);
        for (let i = 0; i < this.dim; i++) {
            newEnd[i] = this.end.getPoint()[i] + that.coords[i];
        }
        return new Vector(this.start, new Point(newEnd), null);
    }

    dot(that) {
        let sum = 0;
        for (let i = 0; i < this.dim; i++) {
            sum += this.coords[i] * that.coords[i];
        }
        return sum;
    }

    len() {
        let sum = 0;
        for (let i = 0; i < this.dim; i++) {
            sum += this.coords[i] ** 2;
        }
        return Math.sqrt(sum);
    }

    norm() {
        let length = this.len();
        let newCoords = [];
        let newEnd = [];
        for (let i = 0; i < this.dim; i++) {
            newCoords[i] = this.coords[i] / length;
            newEnd[i] = this.start.getPoint()[i] + newCoords[i];
        }
        return new Vector(this.start, new Point(newEnd), null);
    }

    scale(scalar) {
        let newCoords = [];
        let newEnd = [];
        for (let i = 0; i < this.dim; i++) {
            newCoords[i] = this.coords[i] * scalar;
            newEnd[i] = this.start.getPoint()[i] + newCoords[i];
        }
        return new Vector(this.start, new Point(newEnd), null);
    }

    perp() {
        let sum = 0;
        for (let i = 0; i < this.dim; i++) {
            sum += this.coords[i];
        }
        let v = - (sum) / this.coords[this.coords.length - 1];
        let coords = []
        for (let i = 0; i < this.dim; i++) {
            coords[i] = 1;
        }
        coords[coords.length - 1] = v;
        let vec = new Vector(null, null, coords).norm();
        return vec;
    }

    getDim() {
        return this.dim;
    }

    travel(point) {
        let end = [];
        for (let i = 0; i < this.dim; i++) {
            end[i] = point.getPoint()[i] + this.coords[i];
        }
        return new Point(end);
    }

    angle(that) {
        return Math.acos(this.dot(that) / (this.len() * that.len()));
    }

    cross(that) {
        let coords1 = this.coords;
        let coords2 = that.coords;
        let x = coords1[1] * coords2[2] - coords1[2] * coords2[1];
        let y = coords1[2] * coords2[0] - coords1[0] * coords2[2];
        let z = coords1[0] * coords2[1] - coords1[1] * coords2[0];
        return new Vector(null, null, [x, y, z]);
    }

    minDist(point) {
        let pl = this.getStart();
        let c = this.dot(new Vector(null, point, null).sub(new Vector(null, pl, null))) / this.dot(this);
        let closest = this.scale(c).travel(pl);
        return point.dist(closest);
    }

    perpVec() {
        let coords = this.coords;
        coords[0] += 1;
        let vec = new Vector(null, null, coords);
        return this.cross(vec);
    }
}

class rotationalMatrix {
    matrix;
    constructor(alpha, beta, gamma) {
        this.matrix = [];
        let m11 = Math.cos(beta) * Math.cos(gamma);
        let m12 = Math.sin(alpha) * Math.sin(beta) * Math.cos(gamma) - Math.cos(alpha) * Math.sin(gamma);
        let m13 = Math.cos(alpha) * Math.sin(beta) * Math.cos(gamma) + Math.sin(alpha) * Math.sin(gamma);

        let m21 = Math.cos(beta) * Math.sin(gamma);
        let m22 = Math.sin(alpha) * Math.sin(beta) * Math.sin(gamma) + Math.cos(alpha) * Math.cos(gamma);
        let m23 = Math.cos(alpha) * Math.sin(beta) * Math.sin(gamma) - Math.sin(alpha) * Math.cos(gamma);

        let m31 = - Math.sin(beta);
        let m32 = Math.sin(alpha) * Math.cos(beta);
        let m33 = Math.cos(alpha) * Math.cos(beta);

        this.matrix[0] = [m11, m12, m13];
        this.matrix[1] = [m21, m22, m23];
        this.matrix[2] = [m31, m32, m33];
    }
    

    vecMul(vector) {
        let start = vector.getStart();
        let end = vector.getEnd();
        let newStart = [];
        let newEnd = [];
        for (let i = 0; i < this.matrix.length; i++) {
            let sum1 = 0;
            let sum2 = 0;
            for (let j = 0; j < vector.getDim(); j++) {
                sum1 += start.getPoint()[j] * this.matrix[i][j];
                sum2 += end.getPoint()[j] * this.matrix[i][j];
            }
            newStart[i] = sum1;
            newEnd[i] = sum2;
        }
        return new Vector(new Point(newStart), new Point(newEnd), null);
    }

    pntMul(point) {
        let newPoint = [];
        for (let i = 0; i < this.matrix.length; i++) {
            let sum = 0;
            for (let j = 0; j < point.getDim(); j++) {
                sum += point.getPoint()[j] * this.matrix[i][j];
            }
            newPoint[i] = sum;
        }
        return new Point(newPoint);
    }
}

class Torus {
    R;
    r;
    center;
    vector;
    matrix;
    spheres;
    constructor(R, r, center, vector) {
        this.R = R;
        this.r = r;
        this.center = center;
        this.vector = vector;
        // this.spheres = [];
        // for (let i = 0; i < spheres; i++) {
        //     let angle = i * 2 * Math.PI / spheres;
        //     let vec = this.vector.perpVec();
        //     this.spheres[i] = new Sphere()
        // }
    }

    rotate(alpha, beta, gamma) {
        this.matrix = new rotationalMatrix(alpha, beta, gamma);
        this.vector = this.matrix.vecMul(this.vector);
    }

    lineCollide(line) {
        let min = line.minDist(this.center);
        if (min > this.R + r) {
            return char;
        }
        let matrix = new rotationalMatrix(rotations * alpha, rotations * beta, rotations * gamma);
        let ch = char;
        let minLen = 200;
        for (let i = 0; i < spheres; i++) {   
            let angle = i * 2 * Math.PI / spheres;
            let vec = new Vector(null, null, [Math.cos(angle), Math.sin(angle), 0]);
            let Rvec = matrix.vecMul(vec).scale(this.R);
            let p = Rvec.travel(this.center);
            min = line.minDist(p);
            if (min <= this.r) {
                let pl = line.getStart();
                let c = line.dot(new Vector(null, p, null).sub(new Vector(null, pl, null))) / line.dot(line);
                let len = line.len() * c - Math.sqrt(this.r ** 2 - min ** 2);
                if (len < minLen) {
                    minLen = len;
                    var index = ((len - 25) / dist) * ascii.length;
                    ch = ascii.charAt(index);
                }
            }
        }   
        return ch;
    }
}

class Sphere {
    center;
    radius;
    constructor(center, radius) {
        this.center = center;
        this.radius = radius;
    }

    contains(point) {
        return this.center.dist(point) <= this.radius;
    }
    
    intersects(line) {
        let minDist = line.minDist(this.center);
        return minDist <= this.radius;
    }

    intersection(line) {
        let linePoint = line.getStart();
        let centerVec = new Vector(null, this.center, null);
        let lineVec = new Vector(null, linePoint, null);
        let c = line.dot(centerVec.sub(lineVec)) / line.dot(line);
        let minDist = line.minDist(this.center);
        let len = line.len() * c - Math.sqrt(this.radius ** 2 - minDist ** 2);
        return line.scale(c).travel(linePoint);
    }

    distance(line) {
        let linePoint = line.getStart();
        let centerVec = new Vector(null, this.center, null);
        let lineVec = new Vector(null, linePoint, null);
        let c = line.dot(centerVec.sub(lineVec)) / line.dot(line);
        let minDist = line.minDist(this.center);
        let len = line.len() * c - Math.sqrt(this.radius ** 2 - minDist ** 2);
        return len;
    }
}

class SolarSystem {
    center;
    bodies;
    constructor(center, bodies) {
        this.center = center;
        this.bodies = bodies;
    }

    rotate(alpha, beta, gamma) {
        for (let i = 0; i < this.bodies.length; i++) {
            let body = this.bodies[i];
            let radius = body.center.dist(this.center);
            let coeff = 1;
            if (radius != 0) {
                coeff = 25 / radius;
            }
            let matrix = new rotationalMatrix(coeff * alpha, coeff * beta, coeff * gamma);
            let center = body.center;
            center = matrix.pntMul(center);
            this.bodies[i].center = center;
        }
    }

    lineCollide(line) {
        for (let i = 0; i < this.bodies.length; i++) {
            let body = this.bodies[i];
            if (body.intersects(line)) {
                let len = body.distance(line);
                var index = ((len - 25) / dist) * ascii.length;
                if (index < 0) {
                    index = 0;
                }
                // console.log(index);
                return ascii.charAt(index);
            }
        }
        return char;
    }
}
 
const ascii = "$@B%8&WM#oahkbdpqwmZO0QLCJUYXzcvunxrjft/\\|()1{}[]?-_+~<>i!lI;:,\"^`'.";
const char = ascii.charAt(ascii.length - 1);

const R = 15;
const r = 5;
const xResolution = 25;
const yResolution = 25;
const xLimit = 25;
const yLimit = 25;
const xRatio = (2 * xLimit) / xResolution;
const yRatio = (2 * xLimit) / yResolution;
let pixel = 10 * (50 / xResolution);

let torusText = document.getElementById("torus-text");
console.log(torusText);
torusText.style.fontSize = pixel + "px";
torusText.style.lineHeight = pixel + "px";

const dist = 75;
const origin = new Point([0, 0, 0]);
const camPnt = new Point([0, 0, 50]);
const camVec = new Vector(origin, camPnt, null);

let rotations = 0;
let torus = new Torus(R, r, new Point([0, 0, 0]), new Vector(null, null, [0, 0, 1]));
const spheres = 30;

let bodies = [
    new Sphere(origin, 3),
    new Sphere(new Point([15, 0, 0]), 2),
    new Sphere(new Point([12, 5, 0]), 2),
];
let solarSystem = new SolarSystem(origin, bodies);

// const alpha = 0;
// const beta = 0;
const gamma = 2 * Math.PI / 89;
const alpha = 2 * Math.PI / 79;
const beta = 2 * Math.PI / 83;

let shape = torus;

function animation() {
    shape.rotate(alpha, beta, gamma);
    let image = "";
    for (let y = -yLimit; y <= yLimit; y += yRatio) {
        for (let x = -xLimit; x <= xLimit; x += xRatio) {
            if (x != -xLimit) {
                image += " ";
            }
            // Get vector then transform according to rotational matrix
            let ray = new Vector(camPnt, new Point([x, y, 0]), null).norm();
            image += shape.lineCollide(ray);
            // console.log(shape.lineCollide(ray)); 
        }
        image += "\n";
    }
    return image;
}

const output = document.getElementById("torus-text");
var intervalId = setInterval(function() {
    let image = animation();
    output.innerHTML = image;
    rotations++;
}, 25);