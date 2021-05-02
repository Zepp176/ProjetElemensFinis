/*
 *  glfem.c - BOV version
 *  Library for EPL1110 : Finite Elements for dummies
 *
 *  Copyright (C) 2021 UCL-EPL : Vincent Legat
 *  All rights reserved.
 *
 *  GLFW  http://www.glfw.org/ (version utilisée 3.3.2)
 *  BOV   https://git.immc.ucl.ac.be/hextreme/NGP/-/tree/master/deps/BOV
 *
 */
 
 
#include "glfem.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>



double   zoom_init;
double   translate_init[2];
double   w_init;
double   h_init;
double   min_colormap;
double   max_colormap;
float    current_text_color[4] = {1.0,0.0,0.0,1.0};
float    current_color[4]      = {0.0,0.0,0.0,1.0};
float    current_line_color[4] = {0.0,0.0,0.0,1.0};
float 	 current_line_width = 0.005;

typedef bov_window_t glfemWindow;
glfemWindow* theCurrentWindow = NULL;



void glfemSetColor(float color[4]) 
{
    for(int i = 0; i < 4; ++i) {
        current_color[i] = color[i]; }
}

void glfemSetTextColor(float color[4]) 
{
    for(int i = 0; i < 4; ++i) {
        current_text_color[i] = color[i]; }
}

void glfemSetLineWidth(float width) 
{
    current_line_width = width;
}

void glfemSetColorLine(float color[4]) 
{
    for(int i = 0; i < 4; ++i) {
        current_line_color[i] = color[i]; }
}

int glfemGetKey(char theKey)
{
    return (glfwGetKey(theCurrentWindow->self,theKey) == GLFW_PRESS) ;
}


static char theCurrentAction = '0';

char glfemGetAction(void)
{
    if (theCurrentAction != '0')  {
          char theLastCurrentAction = theCurrentAction;
          theCurrentAction = '0';
          return theLastCurrentAction;}
    else  return theCurrentAction;
}

//
// Attention : GLFW ne detecte pas le type du clavier et donc : il ne faut utiliser
// que les lettres communes de AZERTY et QWERTY :-)
//
// Il n’y a pas tant de lettres communes (20 lettres). 
// Dans l’ordre d’apparition : e, r, t, u, i, o, p, s, d, f, g, h, j, k, l, x, c, v, b, n
//

static void glfemKeyCallback(GLFWwindow* self,
            int key, int scancode, int action,int mods) 
{
    bov_window_t* window = (bov_window_t*) glfwGetWindowUserPointer(self);
    if (action==GLFW_PRESS || action==GLFW_REPEAT) {
        switch (key) {
          case GLFW_KEY_ESCAPE :
            glfwSetWindowShouldClose(self,GL_TRUE);
            break;
          case GLFW_KEY_H :
            if(window->help_needed==0) window->help_needed = 1;
            else                       window->help_needed = 0;
            break;
          case GLFW_KEY_R :
            window->param.zoom = zoom_init;
            window->param.translate[0] = translate_init[0];
            window->param.translate[1] = translate_init[1];
            break;}}
    if (key == GLFW_KEY_K && action == GLFW_PRESS)  {theCurrentAction = 'K';}
    if (key == GLFW_KEY_E && action == GLFW_PRESS)  {theCurrentAction = 'E';}
    if (key == GLFW_KEY_S && action == GLFW_PRESS)  {theCurrentAction = 'S';}
    if (key == GLFW_KEY_B && action == GLFW_PRESS)  {theCurrentAction = 'B';}
    if (key == GLFW_KEY_R && action == GLFW_PRESS)  {theCurrentAction = 'R';}
    if (key==GLFW_KEY_ESCAPE) glfwSetWindowShouldClose(self,GL_TRUE);
}
	
	
void glfemWindowCreate(const char *windowName,int w, int h,int n,double *x, double *y)
{
    bov_window_t *window = bov_window_new(w,h, windowName);
    theCurrentWindow = window;
    bov_window_set_color(window, (GLfloat[4]){0.9, 0.9, 0.8, 0.0});
    
    //
    // Defining the current viewport and stores it as the reference
    //
  
    double minX  = femMin(x,n);
    double maxX  = femMax(x,n);
    double minY  = femMin(y,n);
    double maxY  = femMax(y,n);
    double sizeX = (maxX-minX)/1.45;
    double meanX = (maxX+minX)/2.0; 
    double sizeY = (maxY-minY)/1.45;
    double meanY = (maxY+minY)/2.0;
    
    double ratio = (GLfloat) h / (GLfloat) w;
    double size = fmax(sizeX,sizeY);
    double left,right,top,bottom;
    if (ratio > 1.0) {
        left = meanX - size;
        right = meanX + size;
        bottom = meanY - size*ratio;
        top = meanY + size*ratio;  }   
    else {
        left = meanX - size/ratio;
        right = meanX + size/ratio;
        bottom = meanY - size;
        top = meanY + size;  }   
    if ((fabs(top-bottom)) >= (fabs(left-right))) {
        window->param.zoom = 1/(0.45*(fabs(top-bottom))); }   
    else {
        window->param.zoom = 1/(0.38*(fabs(left-right))); }   
  
    window->param.translate[0] = -(left + right) / 2.0;
    window->param.translate[1] = -(bottom + top) / 2.0;
  
    zoom_init = window->param.zoom;
    translate_init[0] = window->param.translate[0];
    translate_init[1] = window->param.translate[1];
    w_init = window->size[0];
    h_init = window->size[1];
    
    //
    // Default call back and help message
    //
  
    glfwSetKeyCallback(window->self, glfemKeyCallback);
    glfemWindowSetHelpMessage((const char[]) {
    "   [esc]   Exit\n"
    "    R      Reset zoom and translation\n"
    "    H      Display/hide keyboard shortcuts\n"});
}

void glfemWindowSetHelpMessage(const char *message)
{
    bov_text_delete(theCurrentWindow->help);
    theCurrentWindow->help = bov_text_new((const GLubyte*)(message),GL_STATIC_DRAW);
    bov_text_set_space_type(theCurrentWindow->help, PIXEL_SPACE);
    bov_text_set_fontsize(theCurrentWindow->help, 20.0f); 
    bov_text_set_boldness(theCurrentWindow->help, 0.1f);
    bov_text_set_outline_width(theCurrentWindow->help, 0.5f);
    bov_text_set_color(theCurrentWindow->help,GLFEM_BLACK);  

}

void glfemWindowResetSize(){
		theCurrentWindow->param.zoom = zoom_init;
		theCurrentWindow->param.translate[0] = translate_init[0];
		theCurrentWindow->param.translate[1] = translate_init[1];
}

void glfemWindowUpdate(){
    float w = theCurrentWindow->param.res[0];
    float w2 = theCurrentWindow->size[0];
	  float ratio = w/w2;
	  
    bov_text_set_fontsize(theCurrentWindow->help, ratio*20.0f);
    bov_text_set_pos(theCurrentWindow->help, (GLfloat[2]){-ratio*20.0f, ratio*(theCurrentWindow->size[1] -30.0f)} );
    bov_window_update(theCurrentWindow);
}

void glfemWindowUpdateAndWait(){
    bov_text_set_pos(theCurrentWindow->help, 
            (GLfloat[2]){-20.0f, theCurrentWindow->size[1] - 30.0f} );
    bov_window_update_and_wait_events(theCurrentWindow);
    
}

void glfemWindowFree() 
{
    bov_window_delete(theCurrentWindow);
}

int glfemWindowShouldClose() 
{
    return bov_window_should_close(theCurrentWindow);
}

void glfemDrawMessage(char *message, double pos[2]){

    float w = theCurrentWindow->param.res[0];
    float w2 = theCurrentWindow->size[0];
	  float ratio = w/w2;
  
    bov_text_t* text = bov_text_new((const GLubyte *)message, GL_STATIC_DRAW);
    text->param =  (bov_text_param_t) {
      .fillColor = {current_text_color[0],current_text_color[1],current_text_color[2],current_text_color[3]},
      .outlineColor = {1.0f ,1.0f, 1.0f, 2.0f},
      .pos = {0.0f, 0.0f},
      .shift = {0.0f, 0.0f},
      .fontSize = ratio*20.0f,
      .boldness = 0.0f,
      .outlineWidth = 0.0f,
      .spaceType = PIXEL_SPACE};
    bov_text_set_pos(text, (GLfloat[2]) {pos[0], pos[1]});  
    bov_text_draw(theCurrentWindow, text);
    bov_text_delete(text);
}


void glfemDrawNodes(double *x, double *y, int n) 
{
    GLfloat (*coord)[2] = malloc(sizeof(coord[0])*n);
    for(int i = 0; i < n; ++i){
      coord[i][0] = x[i];
      coord[i][1] = y[i];}
    bov_points_t* points = bov_points_new(coord,n,GL_STATIC_DRAW);
    bov_points_set_color(points,current_color);
    bov_points_set_width(points, current_line_width/zoom_init);
    bov_points_draw(theCurrentWindow,points, 0, BOV_TILL_END);
    bov_points_delete(points);
    free(coord);
}

void glfemDrawElement(double *x, double *y, int n)
{
    GLfloat (*coord)[2] = malloc(sizeof(coord[0])*n);
    for(int i = 0; i < n; ++i){
      coord[i][0] = x[i];
      coord[i][1] = y[i];}
    bov_points_t* points = bov_points_new(coord,n,GL_STATIC_DRAW);
    bov_points_set_color(points,current_color);
    bov_points_set_width(points,current_line_width/zoom_init);
    bov_points_set_outline_width(points, current_line_width/zoom_init);
    bov_points_set_outline_color(points, current_line_color);

    bov_line_loop_draw(theCurrentWindow, points, 0, BOV_TILL_END);
    bov_points_delete(points);
    free(coord);
}

void glfemPlotMesh(femMesh *theMesh)
{
    int i,j,*nodes;
    int nLocalNode = theMesh->nLocalNode;

    GLfloat (*data)[2] = malloc(sizeof(data[0])*nLocalNode*theMesh->nElem);
    
    for (i = 0; i < theMesh->nElem; ++i) {
        nodes = &(theMesh->elem[i*nLocalNode]);
        for (j=0; j < nLocalNode; ++j) {
            data[i*nLocalNode+j][0] = theMesh->X[nodes[j]];
            data[i*nLocalNode+j][1] = theMesh->Y[nodes[j]]; }}
      
    bov_points_t* points = bov_points_new(data,nLocalNode*theMesh->nElem,GL_STATIC_DRAW);
  //   bov_points_set_color(points, GLFEM_BACKGROUND);
    bov_points_set_color(points, current_color);
    bov_points_set_width(points,current_line_width/zoom_init);
    bov_points_set_outline_width(points, 5*current_line_width/zoom_init);
    bov_points_set_outline_color(points, current_line_color);
    
    bov_triangles_draw(theCurrentWindow, points, 0, BOV_TILL_END);
    
    bov_points_delete(points);
    free(data);
}

void glfemSetScale(femMesh *theMesh, double *u)
{
    max_colormap = femMax(u,theMesh->nNode);
    min_colormap = femMin(u,theMesh->nNode);
}

double glfemScale(double minimum, double maximum, double value)
{
    if (value < minimum)        return 0;
    if (minimum == maximum)     return minimum;
    return (value - minimum) / fabs(maximum - minimum);
}

void glfemPlotSolution(femMesh* theMesh, double *u){
    int i,j,*nodes;
    int nLocalNode = theMesh->nLocalNode;
   
    for(int i = 0; i < theMesh->nElem; ++i){
        nodes = &(theMesh->elem[i*nLocalNode]);
        GLfloat (*data)[3] = malloc(sizeof(data[0])*3);
        for (j=0; j < 3; ++j) {
           data[j][0] = theMesh->X[nodes[j]];
           data[j][1] = theMesh->Y[nodes[j]]; 
           data[j][2] = glfemScale(min_colormap,max_colormap,u[nodes[j]]);} 
        bov_points_t* points = bov_points_new_with_value(data,3,GL_STATIC_DRAW);   
        bov_points_set_color(points, GLFEM_BLACK);
        bov_points_set_width(points,current_line_width/zoom_init);
        bov_points_set_outline_color(points, GLFEM_BLACK);
        bov_points_set_outline_width(points, 5*current_line_width/zoom_init);
        bov_triangles_draw(theCurrentWindow, points, 0, BOV_TILL_END);
        bov_points_delete(points);
        free(data);}
   
}


