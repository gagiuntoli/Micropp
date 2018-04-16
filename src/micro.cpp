#include <vector>
#include <boost/program_options.hpp>
#include <iostream>

using namespace boost::program_options;

void on_age(int age)
{
  std::cout << "On age: " << age << '\n';
}

class Problem {

  public:

    int npe = 4;
    int dim = 2;
    int nvoi= 3;

    std::vector<int> bc_nods, bc_y0, bc_y1, bc_x0, bc_x1;
    double xg[4][2];
    double wg[4];
    double b_mat[3][8][4];
    double dsh[4][2][4];
    int nx, ny, nz, nn; 
    double lx, ly, lz, dx, dy, dz;
    int nelem;
    int X0Y0_nod;
    int X1Y0_nod;
    int X1Y1_nod;
    int X0Y1_nod;
    int size_tot;
    int elem_type;
    std::vector<double> strain;
    std::vector<double> stress;
    std::vector<std::vector<int> > elements;
    std::vector<double> int_vars;
    double Ef, Em;
    double Sy_f, Sy_m;

    Problem (int argc, char *argv[]);
};

Problem::Problem (int argc, char *argv[])
{

  try
  {
    options_description desc{"Options Guido"};
    desc.add_options()
      ("help,h", "Help screen")
      ("nx", value<int>()->default_value(10), "Num of Nodes in X dir")
      ("ny", value<int>()->default_value(10), "Num of Nodes in Y dir")
      ("nz", value<int>()->default_value(1) , "Num of Nodes in Z dir")
      ("dim", value<int>()->default_value(2), "Dimension");

    variables_map vm;
    store(parse_command_line(argc, argv, desc), vm);
    notify(vm);

    if (vm.count("help"))
      std::cout << desc << '\n';

    nx = vm["nx"].as<int>();
    ny = vm["ny"].as<int>();
    nz = vm["nz"].as<int>();
    dim = vm["dim"].as<int>();

  }
  catch (const error &ex)
  {
    std::cerr << ex.what() << '\n';
    throw 1;
  }

  lx = 1.0;
  ly = 1.0;
  lz = 1.0;
  dx = lx/(nx-1);
  dy = ly/(ny-1);
  nn = nx*ny*nz;
  if (dim == 2) {
    nelem = (nx-1) * (ny-1);
  }

  wg[0] = 0.25*(dx*dy);
  wg[1] = 0.25*(dx*dy);
  wg[2] = 0.25*(dx*dy);
  wg[3] = 0.25*(dx*dy);

  double xg[4][2] = {
    {-0.577350269189626, -0.577350269189626},
    {+0.577350269189626, -0.577350269189626},
    {+0.577350269189626, +0.577350269189626},
    {-0.577350269189626, +0.577350269189626}};

  for (int gp=0; gp<4; gp++) {
    dsh[0][0][gp] = -(1-xg[gp][1])/4*2/dx ; dsh[0][1][gp] = -(1-xg[gp][0])/4*2/dy;
    dsh[1][0][gp] = +(1-xg[gp][1])/4*2/dx ; dsh[1][1][gp] = -(1+xg[gp][0])/4*2/dy;
    dsh[2][0][gp] = +(1+xg[gp][1])/4*2/dx ; dsh[2][1][gp] = +(1+xg[gp][0])/4*2/dy;
    dsh[3][0][gp] = -(1+xg[gp][1])/4*2/dx ; dsh[3][1][gp] = +(1-xg[gp][0])/4*2/dy;
  }

  for (int gp=0; gp<4; gp++) {
    for (int i=0; i<4; i++) {
      b_mat[1][i*dim][gp] = dsh[i][0][gp]; b_mat[1][i*dim+1][gp] = 0            ;
      b_mat[2][i*dim][gp] = 0            ; b_mat[2][i*dim+1][gp] = dsh[i][1][gp];
      b_mat[3][i*dim][gp] = dsh[i][1][gp]; b_mat[3][i*dim+1][gp] = dsh[i][0][gp];
    }
  }
}
