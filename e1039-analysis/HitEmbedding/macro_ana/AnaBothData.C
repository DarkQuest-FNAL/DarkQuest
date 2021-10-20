R__LOAD_LIBRARY(libana_embedding)
using namespace std;

void AnaBothData()
{
  AnaCleanAndMessyData* ana = new AnaCleanAndMessyData();
  //ana->Verbosity(10);
  ana->Init("not_embedded/ana_tree.root", "embedded/ana_tree.root");
  ana->Analyze();
  ana->End();
  exit(0);
}
