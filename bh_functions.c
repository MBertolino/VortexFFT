#include <stdlib.h>
#include <math.h>
#include "fft_functions.h"
#include "bh_functions.h"
#include "misc.h"

/* Insert a new node */
void insert(node_t **node, double origo_x, double origo_y, double width, double vx, double vy, double x, double y, double dx, double dy, int index) {

  // If null, this is a leaf node
  if ((*node) == NULL) {
    *node = (node_t*)malloc(sizeof(node_t));
    (*node)->origo_x = origo_x;
    (*node)->origo_y = origo_y;
    (*node)->width = width;
    (*node)->vx = vx;
    (*node)->vy = vy;
    (*node)->x = x;
    (*node)->y = y;
    (*node)->dx = dx;
    (*node)->dy = dy;
    (*node)->size = 1;
    (*node)->index = index;
    (*node)->nw = NULL;
    (*node)->ne = NULL;
    (*node)->sw = NULL;
    (*node)->se = NULL;

  // If this is a leaf node, split the node
  } else if ((*node)->size == 1) {
	
    // Add previous particle
    if ((*node)->x < (*node)->origo_x) {
      if ((*node)->y < (*node)->origo_y) // South West
        insert(&((*node)->sw), origo_x - width/4, origo_y - width/4, width/2, (*node)->vx, (*node)->vy,\
			(*node)->x, (*node)->y, (*node)->dx, (*node)->dy, (*node)->index);
      else // North West
        insert(&((*node)->nw), origo_x - width/4, origo_y + width/4, width/2, (*node)->vx, (*node)->vy,\
			(*node)->x, (*node)->y, (*node)->dx, (*node)->dy, (*node)->index);
    } else {
      if ((*node)->y < (*node)->origo_y) // South East
        insert(&((*node)->se), origo_x + width/4, origo_y - width/4, width/2, (*node)->vx, (*node)->vy,\
			(*node)->x, (*node)->y, (*node)->dx, (*node)->dy, (*node)->index);
      else // North East
        insert(&((*node)->ne), origo_x + width/4, origo_y + width/4, width/2, (*node)->vx, (*node)->vy,\
			(*node)->x, (*node)->y, (*node)->dx, (*node)->dy, (*node)->index);
    }
    
    // Add new particle
    if (x < (*node)->origo_x) {
      if (y < (*node)->origo_y) // South West
        insert(&((*node)->sw), origo_x - width/4, origo_y - width/4, width/2, (*node)->vx, (*node)->vy,\
			x, y, dx, dy, index);
      else // North West
        insert(&((*node)->nw), origo_x - width/4, origo_y + width/4, width/2, (*node)->vx, (*node)->vy,\
			x, y, dx, dy, index);
    } else {
      if (y < (*node)->origo_y) // South East
        insert(&((*node)->se), origo_x + width/4, origo_y - width/4, width/2, (*node)->vx, (*node)->vy,\
		       	x, y, dx, dy, index);
      else // North East
        insert(&((*node)->ne), origo_x + width/4, origo_y + width/4, width/2, (*node)->vx, (*node)->vy,\
			x, y, dx, dy, index);
    }
    
    double abs_dx, abs_dy, abs_ndx, abs_ndy;
    abs_dx = fabs(dx);
    abs_dy = fabs(dy);
    abs_ndx = fabs((*node)->dx);
    abs_ndy = fabs((*node)->dy);

    // Update this nodes properties
    (*node)->size += 1;
    //(*node)->vx = ((*node)->vx + vx)*0.5;
    //(*node)->vy = ((*node)->vy + vy)*0.5;
    (*node)->x = ((*node)->x*abs_ndx + x*abs_dx)*0.5/(abs_ndx + abs_dx);
    (*node)->y = ((*node)->y*abs_ndy + y*abs_dy)*0.5/(abs_ndy + abs_dy);
    (*node)->dx += dx;
    (*node)->dy += dy;
		
    // If this is a cluster node
  } else {
    
    // Add new particle
    if (x < (*node)->origo_x) {
      if (y < (*node)->origo_y) // South West
        insert(&((*node)->sw), origo_x - width/4, origo_y - width/4, width/2, (*node)->vx, (*node)->vy,\
	  x, y, dx, dy, index);
      else // North West
        insert(&((*node)->nw), origo_x - width/4, origo_y + width/4, width/2, (*node)->vx, (*node)->vy,\
	  x, y, dx, dy, index);
    } else {
      if (y < (*node)->origo_y) // South East
        insert(&((*node)->se), origo_x + width/4, origo_y - width/4, width/2, (*node)->vx, (*node)->vy,\
	  x, y, dx, dy, index);
      else // North East
        insert(&((*node)->ne), origo_x + width/4, origo_y + width/4, width/2, (*node)->vx, (*node)->vy,\
	  x, y, dx, dy, index);
    }
    
    double abs_dx, abs_dy, abs_ndx, abs_ndy;
    abs_dx = fabs(dx);
    abs_dy = fabs(dy);
    abs_ndx = fabs((*node)->dx);
    abs_ndy = fabs((*node)->dy);

    // Update this nodes properties
    (*node)->size += 1;
    //(*node)->vx = ((*node)->vx + vx)*0.5;
    //(*node)->vy = ((*node)->vy + vy)*0.5;
    (*node)->x = ((*node)->x*abs_ndx + x*abs_dx)*0.5/(abs_ndx + abs_dx);
    (*node)->y = ((*node)->y*abs_ndy + y*abs_dy)*0.5/(abs_ndy + abs_dy);
    (*node)->dx += dx;
    (*node)->dy += dy;
  }
} //*/


/* Barnes-Hut approach */
void vfield_bh(node_t *root, node_t *tree, double theta, double Q, int N, double alpha) {
  if (tree == NULL) return;
  
  if (tree->size != 1) {
    vfield_bh(root, tree->nw, theta, Q, N, alpha);
    vfield_bh(root, tree->ne, theta, Q, N, alpha);
    vfield_bh(root, tree->sw, theta, Q, N, alpha);
    vfield_bh(root, tree->se, theta, Q, N, alpha);
  } else {
    double force_x = 0;  
    double force_y = 0;  
    force_function(root, root, &force_x, &force_y, tree->x, tree->y, tree->dx, tree->dy, tree->index, Q, alpha);
    tree->vx = force_x*theta/(double)N;
    tree->vy = force_y*theta/(double)N;
  
  }
} //*/



/* Compute the total force on a particle */
void force_function(node_t *root, node_t *tree, double *force_x, double *force_y, double x, double y, double dx, double dy,\
		int index, double Q, double alpha) {
  
  if (tree == NULL) return;
  
  // First check if this is a leaf node
  if (tree->size != 1) {
    double q = (tree->width)/distance(x, y, tree->x, tree->y);
    if (q > Q) {
      // Traverse down the tree if the theta condition is not met
      force_function(root, tree->nw, force_x, force_y, x, y, dx, dy, index, Q, alpha);
      force_function(root, tree->ne, force_x, force_y, x, y, dx, dy, index, Q, alpha);
      force_function(root, tree->sw, force_x, force_y, x, y, dx, dy, index, Q, alpha);
      force_function(root, tree->se, force_x, force_y, x, y, dx, dy, index, Q, alpha);
      
      return;
    }
  }
  
  if (tree->index == index)
    return;

  // Compute the relative vector and the distance
  double alpha_d_x = pow(distance(x, y, tree->x, tree->y), alpha);
  //printf("force[%d, %d] = %lf\t%lf\n", index, tree->index, ((tree->dx)-dx)/alpha_d_x, ((tree->dy)-dy)/alpha_d_x); 
  
  // Update the forces
  *force_x += ((tree->dx) - dx)/alpha_d_x;
  *force_y += ((tree->dy) - dy)/alpha_d_x;

  return;
} //*/


/* Recursively write the tree into a file */
void write_tree(node_t *tree, double **dxdt) {
  if (tree == NULL) {
    return;
  } else if (tree->size == 1) {
    (*dxdt)[2*tree->index] = tree->vx;
    (*dxdt)[2*tree->index + 1] = tree->vy;
  } else {
    write_tree(tree->nw, dxdt);
    write_tree(tree->ne, dxdt);
    write_tree(tree->sw, dxdt);
    write_tree(tree->se, dxdt);
  }
}

