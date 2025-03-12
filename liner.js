class Point {
    point;
    constructor(point) {
        this.point = point;
    }

    getPoint() {
        return this.point;
    }

    sub(that) {
        return new Vector(new Point(that.point), new Point(this.point), null);
    }

    dim() {
        return this.point.length;
    }
}

class Vector {
    start;
    end;
    coords;
    // let ray = new Vector(camPnt.getEnd(), new Point([x, y, 0]), null);

    constructor(start, end, coords) {
        if (start != null) {
            this.start = start;
            this.end = end;
            coords = [];
            for (let i = 0; i < this.start.dim(); i++) {
                coords[i] = end.getPoint()[i] - start.getPoint()[i];
            }
            this.coords = coords;
            return;
        }
        start = [];
        for (let i = 0; i < coords.length; i++) {
            start[i] = 0;
        }
        this.start = new Point(start);
        this.end = new Point(coords);
        this.coords = coords;
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
        for (let i = 0; i < this.end.dim(); i++) {
            newEnd[i] = this.end.getPoint()[i] + that.coords[i];
        }
        return new Vector(this.start, new Point(newEnd), null);
    }

    sub(that) {
        let newEnd = [];
        that = that.scale(-1);
        for (let i = 0; i < this.end.dim(); i++) {
            newEnd[i] = this.end.getPoint()[i] + that.coords[i];
        }
        return new Vector(this.start, new Point(newEnd), null);
    }

    dot(that) {
        let sum = 0;
        for (let i = 0; i < this.coords.length; i++) {
            sum += this.coords[i] * that.coords[i];
        }
        return sum;
    }

    len() {
        let sum = 0;
        for (let i = 0; i < this.coords.length; i++) {
            sum += this.coords[i] ** 2;
        }
        return Math.sqrt(sum);
    }

    norm() {
        let length = this.len();
        let newCoords = [];
        let newEnd = [];
        for (let i = 0; i < this.coords.length; i++) {
            newCoords[i] = this.coords[i] / length;
            newEnd[i] = this.start.getPoint()[i] + newCoords[i];
        }
        // console.log(this.coords);
        return new Vector(this.start, new Point(newEnd), null);
    }

    scale(scalar) {
        let newCoords = [];
        let newEnd = [];
        for (let i = 0; i < this.coords.length; i++) {
            newCoords[i] = this.coords[i] * scalar;
            newEnd[i] = this.start.getPoint()[i] + newCoords[i];
        }
        return new Vector(this.start, new Point(newEnd), null);
    }

    perp() {
        let sum = 0;
        for (let i = 0; i < this.coords.length; i++) {
            sum += this.coords[i];
        }
        let v = - (sum) / this.coords[this.coords.length - 1];
        let coords = []
        for (let i = 0; i < this.coords.length; i++) {
            coords[i] = 1;
        }
        coords[coords.length - 1] = v;
        let vec = new Vector(null, null, coords).norm();
        return vec;
    }

    dim() {
        return this.coords.length;
    }
}

class Matrix {
    matrix;
    constructor(matrix, width) {
        if (matrix == null) {
            this.matrix = [];
            for (let i = 0; i < width; i++) {
                this.matrix[i] = [];
                for (let j = 0; j < width; j++) {
                    this.matrix[i][j] = 0;
                }
            }
            return;
        }
        this.matrix = [];
        for (let i = 0; i < matrix.length; i++) {
            this.matrix[i] = [];
            for (let j = 0; j < width; j++) {
                this.matrix[i][j] = matrix[i][j];
            }
        }
    }

    vecMul(vector) {
        let start = vector.getStart();
        let end = vector.getEnd();
        let newStart = [];
        let newEnd = [];
        for (let i = 0; i < this.matrix.length; i++) {
            let sum1 = 0;
            let sum2 = 0;
            for (let j = 0; j < vector.dim(); j++) {
                sum1 += start.getPoint()[j] * matrix[i][j];
                sum2 += end.getPoint()[j] * matrix[i][j];
            }
            newStart[i] = sum1;
            newEnd[i] = sum2;
        }
        return new Vector(new Point(newStart), new Point(newEnd), null);
    }

    pntMul(point) {

    }
}

class rotationalMatrix extends Matrix {
    matrix;
    constructor(alpha, beta, gamma) {
        super(null, 3);
        this.matrix = [];
        let m11 = Math.cos(beta) * Math.cos(gamma);
        let m12 = Math.sin(alpha) * Math.sin(beta) * Math.cos(gamma) - Math.cos(alpha) * Math.sin(gamma);
        let m13 = Math.cos(alpha) * Math.sin(beta) * Math.cos(gamma) + Math.cos(alpha) * Math.sin(gamma);

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
    

    mul(vector) {
        let start = vector.getStart();
        let end = vector.getEnd();
        let newStart = [];
        let newEnd = [];
        for (let i = 0; i < this.matrix.length; i++) {
            let sum1 = 0;
            let sum2 = 0;
            for (let j = 0; j < vector.dim(); j++) {
                sum1 += start.getPoint()[j] * this.matrix[i][j];
                sum2 += end.getPoint()[j] * this.matrix[i][j];
            }
            newStart[i] = sum1;
            newEnd[i] = sum2;
        }
        return new Vector(new Point(newStart), new Point(newEnd), null);
    }
}
 
const ascii = ".'`^\",:;Il!i><~+_-?][}{1)(|\\/tfjrxnuvczXYUJCLQ0OZmwqpdbkhao*#MW&8%B@$";
const char = ascii.charAt(10);

let R = 10;
let r = 5;
let c2 = 2 * Math.PI;
let resolution = 50;
const dist = 50;

const origin = new Point([0, 0, 0]);
const camPnt = new Point([0, 0, 50]);

let rotations = 0;

function brightness(point) {
    var vec = camPnt.sub(point);
    var len = vec.len();
    var index = (len / dist) * ascii.length;
    return ascii.charAt(index);
}

function torus(y) {
    let image = "";
    let rotMatrix = new rotationalMatrix(rotations * Math.PI / 32, rotations * Math.PI / 32, 0);
    let torusVec = rotMatrix.mul(new Vector(null, null, [0, 0, 1]));
    let camVec = rotMatrix.mul(new Vector(null, null, [0, 0, 20]));

    loop:
    for (let x = -resolution; x <= resolution; x++) {
        // Get vector then transform according to rotational matrix
        let ray = new Vector(camVec.getEnd(), new Point([x, y, 0]), null);

        // Project said vector onto the plane defined by the normal vector of the torus
        let proj = ray.sub(torusVec.scale(ray.dot(torusVec)));
        let pos = proj.norm().scale(R);
        let neg = proj.norm().scale(-R);
        if (ray.sub(neg).len() < ray.sub(pos).len()) {
            proj = neg;
        }
        else {
            proj = pos;
        }

        // ray = ray.norm().scale(100);
        for (let i = 0; i < 150; i++) {
            ray = ray.sub(ray.norm().scale(r / 15));
            let diff = ray.sub(proj);
            if (diff.len() < r) {
                let end = ray.getEnd();
                while (diff.len() < r) {
                    ray = ray.sub(ray.norm().scale(r / 30));
                    end = ray.getEnd();
                }
                let c = brightness(end);
                image += c;
                continue loop;
            }
        }
        image += char;
        if (x != resolution - 1) {
            image += " ";
        }
    }
    image += "\n";
    return image;
}

onmessage = (e) => {
    let y = e.data[0];
    let line = torus(y);
    rotations++;
    postMessage(line);
}