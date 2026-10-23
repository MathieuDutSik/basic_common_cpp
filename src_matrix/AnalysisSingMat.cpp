// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>

// clang-format off
#include "NumberTheory.h"
#include "MAT_Matrix.h"
#include "MAT_MatrixInt.h"
#include "Namelist.h"
// clang-format on

FullNamelist NAMELIST_GetStandardAnalysis() {
  std::map<std::string, SingleBlock> ListBlock;
  // PROC
  {
    std::map<std::string, std::string> ListIntValues_doc;
    std::map<std::string, std::string> ListBoolValues_doc;
    std::map<std::string, std::string> ListDoubleValues_doc;
    std::map<std::string, std::string> ListStringValues_doc;
    ListStringValues_doc["InputFile"] =
        "The file containing the real symmetric matrix M";
    ListDoubleValues_doc["thrEigenvalueZero"] = "Default: 0.0001\n\
The eigenvalues of smaller absolute value are considered to be zero";
    ListIntValues_doc["MaximumValueIntSearch"] = "Default: 500\n\
The bound of the search for the integral approximations of the vectors";
    ListBoolValues_doc["SearchKernel"] = "Default: F\n\
Whether to search for an integral basis of the kernel";
    ListBoolValues_doc["SearchOrthKernel"] = "Default: F\n\
Whether to search for an integral basis of the orthogonal of the kernel";
    ListBoolValues_doc["ShowKernel"] = "Default: F\n\
Whether to print the eigenvectors of the kernel";
    ListBoolValues_doc["ShowOrthKernel"] = "Default: F\n\
Whether to print the eigenvectors of the orthogonal of the kernel";
    ListBoolValues_doc["CanonicalizeByMinValue"] = "Default: F\n\
Whether to divide the printed vectors by their smallest entry in absolute\n\
value";
    ListBoolValues_doc["CheckPrimalDualCancellation"] = "Default: F\n\
Whether to print the trace of M P, P the matrix of PrimalDualPairFile";
    ListStringValues_doc["PrimalDualPairFile"] = "Default: unset\n\
The file of the matrix P, used if CheckPrimalDualCancellation = T";
    ListBlock["PROC"].setListIntValues_doc(ListIntValues_doc);
    ListBlock["PROC"].setListBoolValues_doc(ListBoolValues_doc);
    ListBlock["PROC"].setListDoubleValues_doc(ListDoubleValues_doc);
    ListBlock["PROC"].setListStringValues_doc(ListStringValues_doc);
  }
  return FullNamelist(ListBlock);
}

int main(int argc, char *argv[]) {
  HumanTime time;
  std::cerr << std::setprecision(20);
  try {
    FullNamelist eFull = NAMELIST_GetStandardAnalysis();
    if (argc != 2) {
      std::cerr << "This program is used as\n";
      std::cerr << "AnalysisSingMat [Namelist]\n";
      eFull.NAMELIST_WriteNamelistFile(std::cerr, true);
      return -1;
    }
    std::string eFileName = argv[1];
    NAMELIST_ReadNamelistFile(eFileName, eFull);
    SingleBlock const &BlPROC = eFull.get_block("PROC");
    std::string InputFile = BlPROC.get_string("InputFile");
    double thrEigenvalueZero = BlPROC.get_double("thrEigenvalueZero");
    double MaximumValueIntSearch = BlPROC.get_int("MaximumValueIntSearch");
    bool SearchKernel = BlPROC.get_bool("SearchKernel");
    bool SearchOrthKernel = BlPROC.get_bool("SearchOrthKernel");
    bool ShowOrthKernel = BlPROC.get_bool("ShowOrthKernel");
    bool ShowKernel = BlPROC.get_bool("ShowKernel");
    bool CanonicalizeByMinValue = BlPROC.get_bool("CanonicalizeByMinValue");
    bool CheckPrimalDualCancellation =
        BlPROC.get_bool("CheckPrimalDualCancellation");
    std::string PrimalDualPairFile = BlPROC.get_string("PrimalDualPairFile");
    // reading the matrix
    if (!FILE_IsExistingFile(InputFile)) {
      std::cerr << "The InputFile=" << InputFile << " is missing\n";
      throw TerminalException{1};
    }
    std::ifstream INmat(InputFile);
    MyMatrix<double> eMat = ReadMatrix<double>(INmat);
    int nbRow = eMat.rows();
    if (CheckPrimalDualCancellation) {
      if (!FILE_IsExistingFile(PrimalDualPairFile)) {
        std::cerr << "The PrimalDualPairFile=" << PrimalDualPairFile
                  << " is missing\n";
        throw TerminalException{1};
      }
      std::ifstream PAIRmat(PrimalDualPairFile);
      std::cerr << "PrimalDualPairFile=" << PrimalDualPairFile << "\n";
      MyMatrix<double> PairMat = ReadMatrix<double>(PAIRmat);
      MyMatrix<double> eProd = eMat * PairMat;
      double eScal = eProd.trace();
      std::cerr << "PairDual signature = " << eScal << "\n";
    }

    std::cerr << "Read eMat rows/cols = " << eMat.rows() << " / " << eMat.cols()
              << "\n";
    //
    double sumCoeff = 0;
    for (int iRow = 0; iRow < nbRow; iRow++) {
      for (int iCol = 0; iCol < nbRow; iCol++) {
        sumCoeff += fabs(eMat(iRow, iCol));
      }
    }
    std::cerr << "sumCoeff = " << sumCoeff << "\n";
    // creating the RecSparse operators.
    Eigen::SelfAdjointEigenSolver<MyMatrix<double>> eig(eMat);
    MyVector<double> ListEig = eig.eigenvalues();
    MyMatrix<double> ListVect = eig.eigenvectors();
    std::cerr << "ListEigenvalues =";
    for (int iRow = 0; iRow < nbRow; iRow++)
      std::cerr << " " << ListEig(iRow);
    std::cerr << "\n";
    int DimKernel = 0;
    std::vector<MyVector<double>> ListKernel;
    std::vector<MyVector<double>> ListOrthKernel;
    WriteVector(std::cerr, ListEig);
    auto TheCan = [&](MyVector<double> const &V) -> MyVector<double> {
      if (!CanonicalizeByMinValue)
        return V;
      double eVal = V.cwiseAbs().minCoeff();
      std::cerr << "eVal=" << eVal << "\n";
      return V / eVal;
    };
    for (int iRow = 0; iRow < nbRow; iRow++) {
      double eEig = ListEig(iRow);
      MyVector<double> eCol = GetMatrixCol(ListVect, iRow);
      MyVector<double> TheDiff = eEig * eCol - eMat * eCol;
      //      double err=TheDiff.cwiseAbs().maxCoeff();
      //      std::cerr << "iRow=" << iRow << " eEig=" << eEig << " err=" << err
      //      << "\n";
      if (fabs(eEig) < thrEigenvalueZero) {
        DimKernel++;
        ListKernel.push_back(eCol);
        if (ShowKernel)
          WriteVector(std::cerr, TheCan(eCol));
      } else {
        ListOrthKernel.push_back(eCol);
        if (ShowOrthKernel) {
          WriteVector(std::cerr, TheCan(eCol));
        }
      }
    }
    using Tint = mpz_class;
    auto Reduction = [&](std::vector<MyVector<double>> const &inpListVect)
        -> std::pair<double, MyMatrix<Tint>> {
      MyMatrix<double> inpMat = MatrixFromVectorFamily(inpListVect);
      MyMatrix<double> RedMat = CanonicalizeBasisVectorSpace(inpMat);
      std::vector<MyVector<Tint>> LVect;
      int nbRow = RedMat.rows();
      double SumErr = 0;
      for (int iRow = 0; iRow < nbRow; iRow++) {
        MyVector<double> eRow = GetMatrixRow(RedMat, iRow);
        std::pair<double, MyVector<Tint>> PairRed =
            FindBestIntegerApproximation<Tint, double>(eRow,
                                                       MaximumValueIntSearch);
        LVect.push_back(PairRed.second);
        SumErr += PairRed.first;
      }
      return {SumErr, MatrixFromVectorFamily(LVect)};
    };
    if (SearchKernel) {
      std::pair<double, MyMatrix<Tint>> PairRed = Reduction(ListKernel);
      std::cerr << "Approximate kernel found with error=" << PairRed.first
                << "\n";
      WriteMatrix(std::cerr, PairRed.second);
    }
    if (SearchOrthKernel) {
      std::pair<double, MyMatrix<Tint>> PairRed = Reduction(ListOrthKernel);
      std::cerr << "Approximate orthogonal of kernel found with error="
                << PairRed.first << "\n";
      WriteMatrix(std::cerr, PairRed.second);
    }
    std::cerr << "DimKernel=" << DimKernel << " nbRow=" << nbRow << "\n";
  } catch (TerminalException const &e) {
    exit(e.eVal);
  }
  runtime(time);
}
