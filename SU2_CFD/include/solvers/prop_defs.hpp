#define BEM_MAXR      50
#define BEM_MAXALF   100
#define BEM_MAX_ITER  20

typedef struct 
{
  int nblades; 
  float dia; 
  float rhub; 
  float ang0_75; 
} propeller_geom_struct; 

typedef struct 
{
  int nblades; 
  su2double dia; 
  su2double rhub; 
  su2double ang0_75; 
} dpropeller_geom_struct; 

typedef struct 
{
  int nalf; 
  int nrad; 
  float r1[BEM_MAXR]; 
  float chord[BEM_MAXR]; 
  float setangle[BEM_MAXR]; 
  float alf[BEM_MAXALF][BEM_MAXR]; 
  float cl_arr[BEM_MAXALF][BEM_MAXR]; 
  float cd_arr[BEM_MAXALF][BEM_MAXR]; 
} propeller_section_struct; 

typedef struct 
{
  int nalf; 
  int nrad; 
  su2double r1[BEM_MAXR]; 
  su2double chord[BEM_MAXR]; 
  su2double setangle[BEM_MAXR]; 
  su2double alf[BEM_MAXALF][BEM_MAXR]; 
  su2double cl_arr[BEM_MAXALF][BEM_MAXR]; 
  su2double cd_arr[BEM_MAXALF][BEM_MAXR]; 
} dpropeller_section_struct; 
