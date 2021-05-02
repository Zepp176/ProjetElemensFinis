static GLchar points_vert[]={"#version 150 core\n"
"in vec2 pos;\n"
"in float value;\n"
"out float vertexvalue;\n"
"void main() { \n"
"  gl_Position = vec4(pos, 0.0, 1.0);\n"
"  vertexvalue = value;\n"
"}\n"
};
