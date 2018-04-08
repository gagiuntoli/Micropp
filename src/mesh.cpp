#include <iostream>
#include <string>
#include <vector>

using namespace std;

class Node {

  public:

    Node () {};
    Node (vector<double> Coord) {this->Coord = Coord;};
    vector<double> GetCoord () {return this->Coord;};
    int GetNum () {return this->Num;};
    int SetCoord (vector<double> Coord) {this->Coord = Coord;};
    int SetNum (int Num) {this->Num = Num;};

  private:

    int Num;
    vector<double> Coord;

};

class Element {

  public:

    Element (vector<Node> Nodes);
    vector<Node> GetNodes () {return this->Nodes;};
    string GetType () {return this->Type;};

  private:

    int Num;
    int nNods;
    string Type;
    vector<Node> Nodes;

};

Element::Element (vector<Node> Nodes)
{
  this->Nodes = Nodes;
  this->nNods = Nodes.size();
  if (Nodes.size() == 4) {
    this->Type = "QUAD4";
  }
}

class Mesh {

  public:

    Mesh (int nNodsX, int nNodsY);
    Mesh (int nNodsX, int nNodsY, int nNodsZ);
    vector<int> GetNumNods ();
    void Enumerate (string method);

  private:

    int Dimension;
    int nNodsX, nNodsY, nNodsZ;
    int nElemX, nElemY, nElemZ;
    int nNods, nElems;
    double Lx, Ly, Lz;
    double hx, hy, hz;
    vector<Element> Elements;

};

Mesh::Mesh (int nNodsX, int nNodsY)
{ 
  this->nNodsX = nNodsX;
  this->nNodsY = nNodsY;
  this->nElemX = this->nNodsX - 1;
  this->nElemY = this->nNodsY - 1;
  this->Lx = 1.0;
  this->Ly = 1.0;
  this->hx = this->Lx / this->nNodsX;
  this->hy = this->Ly / this->nNodsY;
  this->Dimension = 2;

  for (int ex=0; ex<this->nElemX; ++ex) {
    for (int ey=0; ey<this->nElemY; ++ey) {
      vector<Node> nodes (4);
      // complete nodes vector
      // create one element
      Element elem (nodes);
      this->Elements.push_back(elem); // insert in the vector
    }
  }

}

Mesh::Mesh (int nNodsX, int nNodsY, int nNodsZ)
{ 
  this->nNodsX = nNodsX;
  this->nNodsY = nNodsY;
  this->nNodsZ = nNodsZ;
  this->Dimension = 3;
}

vector<int> Mesh::GetNumNods (void)
{
  vector<int> v;
  if (this->Dimension == 2) {
    v.push_back(nNodsX);
    v.push_back(nNodsY);
  } else if (this->Dimension == 3) {
    v.push_back(nNodsX);
    v.push_back(nNodsY);
    v.push_back(nNodsZ);
  }
  return v;
}

void Mesh::Enumerate (string method)
{
  if (method == "easy") {
    if (this->Dimension == 2) {

    } else if (this->Dimension == 3) {
    }
  }
}

int main (int argc, char *argv[])
{
  Mesh mesh2d (2, 4);
  Mesh mesh3d (5, 4, 3);

  cout << "Testing Mesh Class" << endl;

  vector<int> vec = mesh2d.GetNumNods();
  for (vector<int>::const_iterator i = vec.begin(); i != vec.end(); ++i)
    cout << *i << ' ';
  cout << endl;

  vec = mesh3d.GetNumNods();
  for (vector<int>::const_iterator i = vec.begin(); i != vec.end(); ++i)
    cout << *i << ' ';
  cout << endl;


  cout << "Testing Node Class" << endl;
  vector<double> vec_d {1.2, 4.2, 5.1};
  Node nod_1 (vec_d);
  vector<double> vec_e = nod_1.GetCoord();
  for (vector<double>::const_iterator i = vec_e.begin(); i != vec_e.end(); ++i)
    cout << *i << ' ';
  cout << endl;

  cout << "Testing Element Class" << endl;
  vec_d = {-1.2, 8.2, 5.1};
  Node nod_2 (vec_d);
  vec_d = {1.2, 0.2, 5.1};
  Node nod_3 (vec_d);
  vec_d = {3.2, 4.2001, 54.1};
  Node nod_4 (vec_d);
  vector<Node> nod_vec;
  nod_vec.push_back(nod_1);
  nod_vec.push_back(nod_2);
  nod_vec.push_back(nod_3);
  nod_vec.push_back(nod_4);
  Element elem_1 (nod_vec);
  vector<Node> vec_c = elem_1.GetNodes();

  return 0;
}
