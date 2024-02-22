const canvas = document.getElementById("c");
const gl = canvas.getContext("webgl2");
if (!gl) {
	console.error("WebGL2 is not supported.");
}

const vertexShaderSource = `#version 300 es
	in vec2 a_position;
	uniform vec4 u_color;
	out vec4 v_color;
	void main() {
		gl_Position = vec4(a_position, 0.0, 1.0);
		v_color = u_color;
	}
`;

const fragmentShaderSource = `#version 300 es
	precision highp float;
	in vec4 v_color;
	out vec4 fragColor;
	void main() {
		fragColor = v_color;
	}
`;

const xlim = [-0.5, 1];
const x_offset = 0.07;
const ylim = [-1, 1];
const y_offset = 0.1;

function lmap(v, vlim) {
	// Linear map from vlim[-1] => -1 and vlim[1] => 1
	const a = vlim[0];
	const b = vlim[1];
	return (-2/(a-b))*v + (a/(a-b)) + (b/(a-b));
}

const pointVertexCount = 8;
const pointRadius = 8;
const p = [0, 0];

function addPointVertices(p, r=pointRadius, vertexCount=pointVertexCount) {
	const vertices = [];
	for (let i = 0; i < vertexCount; i++) {
		const angle = (i/vertexCount) * Math.PI * 2;
		const x = lmap(p[0], xlim) + r*Math.cos(angle)/gl.canvas.width;
		const y = lmap(p[1], ylim) + r*Math.sin(angle)/gl.canvas.height;
		vertices.push(x, y);
	}
	vertices.push(vertices[0], vertices[1]);

	return vertices;
}

let points = [];
const A = [
	2, 1,
	1, 3
];
const b = [1, 2];
const x_true = [0.2, 0.6];
points.push(x_true);
const x0 = [0.5, 0.5];
points.push(x0);

// Gauss-Seidel
let x = structuredClone(x0);
for (let iter = 0; iter < 10; iter++) {
	let i = iter % 2;
	if (i == 0) {
		let xx = (b[0] - A[1]*x[1])/A[0];
		x[0] = xx;
	} else {
		let xx = (b[1] - A[2]*x[0])/A[3];
		x[1] = xx;
	}
	points.push(structuredClone(x));
}

// Determine limits
xlim[0] = Math.min(...points.map((p) => p[0]));
xlim[1] = Math.max(...points.map((p) => p[0]));
ylim[0] = Math.min(...points.map((p) => p[1]));
ylim[1] = Math.max(...points.map((p) => p[1]));
// Offset
dx = (xlim[1]-xlim[0]) * x_offset;
dy = (ylim[1]-ylim[0]) * y_offset;
xlim[0] -= dx;
xlim[1] += dx;
ylim[0] -= dy;
ylim[1] += dy;

// Control
const slider = document.getElementById("iter");
let numPointsDisplayed = 0;
const display = document.getElementById("iter-value");
slider.oninput = function() {
	numPointsDisplayed = parseInt(this.value);
	display.innerHTML = numPointsDisplayed;
	render();
}

function compileShader(gl, source, type) {
	const shader = gl.createShader(type);
	gl.shaderSource(shader, source);
	gl.compileShader(shader);
	if (!gl.getShaderParameter(shader, gl.COMPILE_STATUS)) {
		console.error("Error compiling shader: ", gl.getShaderInfoLog(shader));
		gl.deleteShader(shader);
		return null;
	}
	return shader;
}

const vertexShader = compileShader(gl, vertexShaderSource, gl.VERTEX_SHADER);
const fragmentShader = compileShader(gl, fragmentShaderSource, gl.FRAGMENT_SHADER);

const program = gl.createProgram();
gl.attachShader(program, vertexShader);
gl.attachShader(program, fragmentShader);
gl.linkProgram(program);
if (!gl.getProgramParameter(program, gl.LINK_STATUS)) {
	console.error("Error linking program: ", gl.getProgramInfoLog(program));
}

const vertexBuffer = gl.createBuffer();
gl.bindBuffer(gl.ARRAY_BUFFER, vertexBuffer);

const vao = gl.createVertexArray();
gl.bindVertexArray(vao);
const positionAttribLocation = gl.getAttribLocation(program, "a_position");
gl.vertexAttribPointer(positionAttribLocation, 2, gl.FLOAT, false, 0, 0);
gl.enableVertexAttribArray(positionAttribLocation);
const colorUniformLocation = gl.getUniformLocation(program, "u_color");

function render() {
	gl.clearColor(1.0, 1.0, 1.0, 1.0);
	gl.clear(gl.COLOR_BUFFER_BIT);

	gl.useProgram(program);
	gl.bindVertexArray(vao);

	// Draw lines
	for (let i = 0; i < numPointsDisplayed + 2; i++) {
		// Lines
		gl.uniform4fv(colorUniformLocation, [0.0, 0.0, 0.0, 1.0]);
		if (i > 1) {
			gl.bufferData(gl.ARRAY_BUFFER, new Float32Array([
				lmap(points[i-1][0], xlim), lmap(points[i-1][1], ylim), lmap(points[i][0], xlim), lmap(points[i][1], ylim)
			]), gl.STATIC_DRAW);
			gl.drawArrays(gl.LINE_STRIP, 0, 2);
		}
	}

	// Draw points
	for (let i = 0; i < numPointsDisplayed + 2; i++) {
		if (i == 0) {
			gl.uniform4fv(colorUniformLocation, [0.0, 0.0, 1.0, 1.0]);
		} else {
			gl.uniform4fv(colorUniformLocation, [1.0, 0.0, 0.0, 1.0]);
		}
		gl.bufferData(gl.ARRAY_BUFFER, new Float32Array(addPointVertices(points[i])), gl.STATIC_DRAW);
		gl.drawArrays(gl.TRIANGLE_FAN, 0, pointVertexCount);
	}
}

render();
