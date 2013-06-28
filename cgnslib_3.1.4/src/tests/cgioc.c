/*
   Sample CGIO test program to build files illustrated
   in example database figure.
*/

#include <stdio.h>
#include <ctype.h>
#include <string.h>
#ifdef _WIN32
#include <io.h>
#define unlink _unlink
#else
#include <unistd.h>
#endif

#include "cgns_io.h"

void print_child_list(int cgio_num,double node_id);

int main ()
{
  /* --- Node header character strings */
   char label[CGIO_MAX_LABEL_LENGTH+1];
   char data_type[CGIO_MAX_DATATYPE_LENGTH+1];

  /* --- Database identifier */
  int cgio_num, cgio_num2;

  /* --- Node id variables */
   double root_id,parent_id,child_id,tmp_id,root_id_file2;

  /* --- Data to be stored in database */
   float a[3][4];
   cgsize_t a_dimensions[2] = {4,3};

   int c[6] = {1,2,3,4,5,6};
   cgsize_t c_dimension = 6;

  /* --- miscellaneous variables */
   int i, j;
   int error_state = 1;
   int num_dims, d[6];
   cgsize_t dim_d, dims_b[2];
   float b[3][4];

   for (i = 0; i < 3; i++) {
       for (j = 0; j < 4; j++) {
           a[i][j] = (float)((double)(j + 1) + 0.1 * (double)(i + 1));
       }
   }

  /* ------ begin source code ----- */

  /* --- set database error flag to abort on error */
   cgio_error_abort(error_state);

  /* -------- build file: file_two.cgio ---------- */
  /* --- 1.) open database
         2.) create three nodes at first level
         3.) put label on node f3
         4.) put some data in node f3
         5.) create two nodes below f3
         6.) close database */

   unlink("file_two.cgio");
   cgio_open_file("file_two.cgio",CGIO_MODE_WRITE,CGIO_FILE_NONE,&cgio_num);
   cgio_get_root_id(cgio_num,&root_id);
   root_id_file2 = root_id;
   cgio_create_node(cgio_num,root_id,"f1",&tmp_id);
   cgio_create_node(cgio_num,root_id,"f2",&tmp_id);
   cgio_create_node(cgio_num,root_id,"f3",&parent_id);
   cgio_set_label(cgio_num,parent_id,"label on node f3");

   cgio_set_dimensions(cgio_num,parent_id,"R4",2,a_dimensions);
   cgio_write_all_data(cgio_num,parent_id,a);

   cgio_create_node(cgio_num,parent_id,"f4",&child_id);
   cgio_create_node(cgio_num,parent_id,"f5",&child_id);
   cgio_close_file(cgio_num);

  /* -------- build file: file_one.cgio ---------- */
  /* open database and create three nodes at first level */
   unlink("file_one.cgio");
   cgio_open_file("file_one.cgio",CGIO_MODE_WRITE,CGIO_FILE_NONE,&cgio_num);
   cgio_get_root_id(cgio_num,&root_id);
   cgio_create_node(cgio_num,root_id,"n1",&tmp_id);
   cgio_create_node(cgio_num,root_id,"n2",&tmp_id);
   cgio_create_node(cgio_num,root_id,"n3",&tmp_id);

  /* put three nodes under n1 (two regular and one link) */
   cgio_get_node_id(cgio_num,root_id,"n1",&parent_id);
   cgio_create_node(cgio_num,parent_id,"n4",&tmp_id);
   cgio_create_link(cgio_num,parent_id,"l3","file_two.cgio","/f3",&tmp_id);
   cgio_create_node(cgio_num,parent_id,"n5",&tmp_id);

  /* put two nodes under n4 */
   cgio_get_node_id(cgio_num,parent_id,"n4",&child_id);
   cgio_create_node(cgio_num,child_id,"n6",&tmp_id);
   cgio_create_node(cgio_num,child_id,"n7",&tmp_id);

  /* put one nodes under n6 */
   cgio_get_node_id(cgio_num,root_id,"/n1/n4/n6",&parent_id);
   cgio_create_node(cgio_num,parent_id,"n8",&tmp_id);

  /* put three nodes under n3 */
   cgio_get_node_id(cgio_num,root_id,"n3",&parent_id);
   cgio_create_node(cgio_num,parent_id,"n9",&tmp_id);
   cgio_create_node(cgio_num,parent_id,"n10",&tmp_id);
   cgio_create_node(cgio_num,parent_id,"n11",&tmp_id);

  /* put two nodes under n9 */
   cgio_get_node_id(cgio_num,parent_id,"n9",&child_id);
   cgio_create_node(cgio_num,child_id,"n12",&tmp_id);
   cgio_create_node(cgio_num,child_id,"n13",&tmp_id);

  /* put label and data in n13 */
   cgio_set_label(cgio_num,tmp_id,"Label on Node n13");
   cgio_set_dimensions(cgio_num,tmp_id,"I4",1,&c_dimension);
   cgio_write_all_data(cgio_num,tmp_id,c);

  /* put two nodes under n10 (one normal, one link) */
   cgio_get_node_id(cgio_num,root_id,"/n3/n10",&parent_id);
   cgio_create_link(cgio_num,parent_id,"l1"," ","/n3/n9/n13",&tmp_id);
   cgio_create_node(cgio_num,parent_id,"n14",&tmp_id);

  /* put two nodes under n11 (one normal, one link) */
   cgio_get_node_id(cgio_num,root_id,"/n3/n11",&parent_id);
   cgio_create_link(cgio_num,parent_id,"l2"," ","/n3/n9/n13",&tmp_id);
   cgio_create_node(cgio_num,parent_id,"n15",&tmp_id);

  /* ----------------- finished building file_one.cgio ------------- */

  /* ------------- access and print data --------------- */

  /* access data in node f3 (file_two.cgio) through link l3 */
   cgio_get_node_id(cgio_num,root_id,"/n1/l3",&tmp_id);
   cgio_get_label(cgio_num,tmp_id,label);
   cgio_get_data_type(cgio_num,tmp_id,data_type);
   cgio_get_dimensions(cgio_num,tmp_id,&num_dims,dims_b);
   cgio_read_all_data(cgio_num,tmp_id,b);
   printf (" node f3 through link l3:\n");
   printf ("   label       = %s\n",label);
   printf ("   data_type   = %s\n",data_type);
   printf ("   num of dims = %5d\n",num_dims);
   printf ("   dim vals    = %5d %5d\n",(int)dims_b[0],(int)dims_b[1]);
   printf ("   data:\n");
   for (i=0; i<=3; i++)
     {
       for (j=0; j<=2; j++)
         {
           printf("     %10.2f",b[j][i]);
         };
       printf("\n");
     }

  /* access data in node n13 */
   cgio_get_node_id(cgio_num,root_id,"/n3/n9/n13",&tmp_id);
   cgio_get_label(cgio_num,tmp_id,label);
   cgio_get_data_type(cgio_num,tmp_id,data_type);
   cgio_get_dimensions(cgio_num,tmp_id,&num_dims,&dim_d);
   cgio_read_all_data(cgio_num,tmp_id,d);
   printf (" node n13:\n");
   printf ("   label       = %s\n",label);
   printf ("   data_type   = %s\n",data_type);
   printf ("   num of dims = %5d\n",num_dims);
   printf ("   dim val     = %5d\n",(int)dim_d);
   printf ("   data:\n");
   for (i=0; i<=5; i++)
     {
       printf("     %-4d",d[i]);
     }
   printf("\n\n");

  /* access data in node n13 through l1 */
   cgio_get_node_id(cgio_num,root_id,"/n3/n10/l1",&tmp_id);
   cgio_get_label(cgio_num,tmp_id,label);
   cgio_read_all_data(cgio_num,tmp_id,d);
   printf (" node n13 through l1:\n");
   printf ("   label       = %s\n",label);
   printf ("   data:\n");
   for (i=0; i<=5; i++)
     {
       printf("     %-4d",d[i]);
     }
   printf("\n\n");

  /* access data in node n13 through l2 */
   cgio_get_node_id(cgio_num,root_id,"/n3/n11/l2",&tmp_id);
   cgio_get_label(cgio_num,tmp_id,label);
   cgio_read_all_data(cgio_num,tmp_id,d);
   printf (" node n13 through l2:\n");
   printf ("   label       = %s\n",label);
   printf ("   data:\n");
   for (i=0; i<=5; i++)
     {
       printf("     %-4d",d[i]);
     }
   printf("\n\n");

  /* print list of children under root node */
   print_child_list(cgio_num,root_id);

  /* print list of children under n3 */
   cgio_get_node_id(cgio_num,root_id,"/n3",&tmp_id);
   print_child_list(cgio_num,tmp_id);

  /* re-open file_two and get new root id */
   cgio_open_file("file_two.cgio",CGIO_MODE_READ,CGIO_FILE_NONE,&cgio_num2);
   cgio_get_root_id(cgio_num2,&root_id);
   printf (" Comparison of root id:\n");
   printf ("   file_two.cgio original root id = %g\n",root_id_file2);
   printf ("   file_two.cgio new      root id = %g\n",root_id);

   cgio_close_file(cgio_num);
   cgio_close_file(cgio_num2);
   return 0;
}

void print_child_list(int cgio_num, double node_id)
{

/*
   print table of children given a parent node-id
*/
   char node_name[CGIO_MAX_NAME_LENGTH+1];
   int i, num_children, num_ret;

   cgio_get_name(cgio_num,node_id,node_name);
   cgio_number_children(cgio_num,node_id,&num_children);
   printf ("Parent Node Name = %s\n",node_name);
   printf ("  Number of Children = %2d\n",num_children);
   printf ("  Children Names:\n");
   for (i=1; i<=num_children; i++)
     {
       cgio_children_names(cgio_num,node_id,i,1,CGIO_MAX_NAME_LENGTH+1,
           &num_ret,node_name);
       printf ("     %s\n",node_name);
     }
    printf ("\n");
}

