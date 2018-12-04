#include "Usrcat.hpp"
#include "KeepN.hpp"
#include <iostream>
#include <openbabel/mol.h>
#include <openbabel/obconversion.h>

int main(int argc, char *argv[])
{
  if (argc != 4 && argc != 5)
  {
    std::cerr << "Usage: ./usrcat query.sdf candidates.sdf numtokeep "
                 "[outfile.sdf]\n";
    std::cerr << "If outfile.sdf is ommitted, output is sent to the terminal "
                 "for redirection\n";
    exit(-1);
  }
  KeepNAscending<OpenBabel::OBMol> keepn(std::atoi(argv[3]));

  auto ReadMolFromStream = [](OpenBabel::OBConversion &conv,
                              OpenBabel::OBMol &out) {
    if (conv.Read(&out))
    {
      return true;
    }
    return false;
  };

  OpenBabel::OBMol needle, haystackmember;
  std::ifstream inFileNeedle(argv[1]);
  std::ifstream inFileHaystack(argv[2]);
  OpenBabel::OBConversion obconv;
  obconv.SetInAndOutFormats("sdf", "sdf");
  obconv.SetInStream(&inFileNeedle);
  if (!ReadMolFromStream(obconv, needle))
  {
    std::cerr << "Error reading from query file: " << argv[1] << '\n';
    exit(EXIT_FAILURE);
  };

  uint64_t endFilePos{0};
  inFileHaystack.seekg(0, std::ios::end);
  endFilePos = inFileHaystack.tellg();
  inFileHaystack.seekg(0);
  obconv.SetInStream(&inFileHaystack);

  unsigned int counter = 0;
  auto PrintProgress = [&endFilePos, &counter, &inFileHaystack]() {
    std::cerr << "Progress:\t"
              << (inFileHaystack.tellg() / float(endFilePos) * 100) << "% ("
              << counter / 1000 << "k complete)\r";
  };
  Usrcat usrcat;
  USRCATDESCRIPTORS needleDescriptors = usrcat.GetUsrcatValues(needle);
  for (auto &i : needleDescriptors)
    std::cerr << i << ",";
  std::cerr << "\n";

  while (inFileHaystack.tellg() < endFilePos)
  {
    if (!ReadMolFromStream(obconv, haystackmember))
    {
      std::cerr << "Invalid molecule found immediately before file position "
                   "byte offset: "
                << inFileHaystack.tellg() << "\n";
      continue;
    }
    ++counter;
    keepn.insert(haystackmember,
                 usrcat.score(needleDescriptors,
                              usrcat.GetUsrcatValues(haystackmember)));

    if (counter % 1000 == 0)
    {
      PrintProgress();
    }
  }

  std::cerr << "Progress:\t 100% complete - " << counter
            << " molecules processed\n";
  std::ofstream outFile;

  auto AddLabelToMol = [](OpenBabel::OBMol &mol, const std::string &labelname,
                          const std::string &labelvalue) {
    OpenBabel::OBPairData *label = new OpenBabel::OBPairData;
    label->SetAttribute(labelname);
    label->SetValue(labelvalue);
    mol.SetData(label);
  };

  if (argc == 4)
  {
    obconv.SetOutStream(&std::cout);
  }
  else
  {
    outFile.open(argv[4]);
    obconv.SetOutStream(&outFile);
  }
  for (auto i : keepn.best)
  {
    OpenBabel::OBMol outmol(
        i.first); // keepn.best is const, so we cant add to it
    AddLabelToMol(outmol, std::string("USRCATScore"),
                  std::string(argv[1]) + "," + std::string(argv[2]) + "," +
                      needle.GetTitle() + ", " + +", " + i.first.GetTitle() +
                      ", " + std::to_string(i.second));
    obconv.Write(&outmol);
  }
}
