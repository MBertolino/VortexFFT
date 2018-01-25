
/* Define the quad_node struct */
typedef struct quad_node {
	struct quad_node *nw, *ne, *sw, *se;
	double origo_x;
	double origo_y;
	double width;
	double vx;
	double vy;
	double x;
	double y;
	double dx;
	double dy;
	int size, index;
} node_t;


void insert(node_t **node, double origo_x, double origo_y, double width, double vx, double vy, double x, double y, double dx, double dy, int index);


void vfield_bh(node_t *root, node_t *tree, double theta, double Q, int N, double alpha); 


void force_function(node_t *root, node_t *tree, double *force_x, double *force_y, double x, double y, double dx, double dy,\
		int index, double Q, double alpha);


void write_tree(node_t *tree, double **dxdt);
