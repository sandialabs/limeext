
#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_UnitTestRepository.hpp>
#include <Teuchos_GlobalMPISession.hpp>

int main(int argc,char**argv) {
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);
  bool result = Teuchos::UnitTestRepository::runUnitTestsFromMain(argc, argv);
  std::cout << std::endl ;
  return result;
}



