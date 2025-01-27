#include "LHAPDF/LHAPDF.h"
#include "LHAPDF/Reweighting.h"

using namespace LHAPDF;

void RunPDF() {
  std::cout << "start RunPDF" << std::endl;
  // Initialize the PDF set with the given name
  PDF* pdf = mkPDF("NNPDF31_nnlo_as_0118", 0);

  // Get the xF(x) values for the proton at Q=10 GeV
  std::vector<double> xFs(100, 0);

  std::cout << "start loop" << std::endl;
  
  for (int i = 0; i < xFs.size(); i++) {
    xFs[i] = pdf->xfxQ(21, 0.1, 10.0);
  }

  std::cout << "end loop" << std::endl;

  // Print the values
  for (int i = 0; i < xFs.size(); i++) {
    std::cout << xFs[i] << std::endl;
  }
}

int MAXVALUE = pow(2, 32);

double PDFWeight(int id1, int id2, double x1, double x2, double Q2, const PDF& pdf) {
  // handles overflow
  if (id1 > 100)
    id1 = id1 - MAXVALUE;
  if (id2 > 100)
    id2 = id2 - MAXVALUE;

  double w_pdf = pdf.xfxQ2(id1, x1, Q2) * pdf.xfxQ2(id2, x2, Q2);
  return w_pdf;
}

double PDFReweight(int id1, int id2, double x1, double x2, double Q2, const PDF& pdf1, const PDF& pdf2) {
  // handles overflow
  if (id1 > 100)
    id1 = id1 - MAXVALUE;
  if (id2 > 100)
    id2 = id2 - MAXVALUE;

  double wgt = weightxxQ2(id1, id2, x1, x2, Q2, pdf1, pdf2);
  return wgt;
}

void test() {
  std::cout << "start test" << std::endl;
  const PDF* pdf1 = mkPDF("NNPDF31_nnlo_as_0118", 0);
  const PDF* pdf2 = mkPDF("CT18NNLO_as_0118", 0);
  double wgt = PDFReweight(21, 21, 0.1, 0.1, 30, *pdf1, *pdf2);
  std::cout << "wgt = " << wgt << std::endl;
}