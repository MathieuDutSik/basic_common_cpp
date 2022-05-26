// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_SVG_SVGFUNCTIONS_H_
#define SRC_SVG_SVGFUNCTIONS_H_

#include <algorithm>
#include <string>
#include <vector>

struct coor {
  double x;
  double y;
};

coor operator+(coor const &c1, coor const &c2) {
  return {c1.x + c2.x, c1.y + c2.y};
}

coor operator-(coor const &c1, coor const &c2) {
  return {c1.x - c2.x, c1.y - c2.y};
}

std::ostream &operator<<(std::ostream &os, coor const &obj) {
  os << "(" << obj.x << "," << obj.y << ")";
  return os;
}

double scal(coor const &v1, coor const &v2) {
  return v1.x * v2.x + v1.y * v2.y;
}

double norm(coor const &c) { return sqrt(c.x * c.x + c.y * c.y); }

coor MultScal(coor const &c, double const &eScal) {
  return {eScal * c.x, eScal * c.y};
}

coor IsobarycenterPoint(std::vector<coor> const &ListPt) {
  double eX = 0;
  double eY = 0;
  for (auto &ePt : ListPt) {
    eX += ePt.x;
    eY += ePt.y;
  }
  double eF = 1 / double(ListPt.size());
  return {eX * eF, eY * eF};
}

coor RotateCoor(coor const &c, double const &eAngDeg) {
  double eAngRad = eAngDeg * (M_PI / double(180));
  //  std::cerr << "eAngRad=" << eAngRad << "\n";
  double eCos = cos(eAngRad);
  double eSin = sin(eAngRad);
  double cX = c.x * eCos - c.y * eSin;
  double cY = c.x * eSin + c.y * eCos;
  return {cX, cY};
}

std::vector<coor> GetListTangent(std::vector<coor> const &ListPoint) {
  int len = ListPoint.size();
  std::vector<coor> ListTangent(len);
  for (int i = 0; i < len; i++) {
    coor eTangent;
    if (i == 0) {
      eTangent = ListPoint[1] - ListPoint[0];
    } else {
      if (i == len - 1) {
        eTangent = ListPoint[len - 1] - ListPoint[len - 2];
      } else {
        coor coorSum = ListPoint[i + 1] - ListPoint[i - 1];
        eTangent = MultScal(coorSum, double(1) / double(2));
      }
    }
    ListTangent[i] = eTangent;
  }
  return ListTangent;
}

struct SVGclippath {
  std::string name;
  double x;
  double y;
  double width;
  double height;
};

struct SVGqualInfo {
  std::vector<int> color;
  double Size;
  std::string MarkerEnd;
  std::string clip;
};

struct SVGqualInfoPolyline {
  std::vector<int> colorfill;
  std::vector<int> colorstroke;
  double Size;
  std::string MarkerEnd;
  std::string clip;
};

struct SVGbezier {
  coor pointM;
  coor pointC;
  coor point1;
  coor point2;
  SVGqualInfo eQual;
};

struct SVGtext {
  coor point;
  std::string string;
  std::vector<int> color;
  int Size;
};

struct SVGpolyline {
  std::vector<coor> ListCoor;
  SVGqualInfoPolyline eQual;
};

struct SVGline {
  coor ePt;
  coor fPt;
  SVGqualInfo eQual;
};

struct SVGellipse {
  coor c;
  coor r;
  SVGqualInfo eQual;
};

enum SvgType { clip, polyline, bezier, line, ellipse, text };

struct SVGgeneral {
public:
  SVGgeneral(SVGclippath const &eClip) {
    svgtype = SvgType::clip;
    clip = eClip;
  }
  SVGgeneral(SVGpolyline const &ePolyline) {
    svgtype = SvgType::polyline;
    polyline = ePolyline;
  }
  SVGgeneral(SVGbezier const &eBez) {
    svgtype = SvgType::bezier;
    bezier = eBez;
  }
  SVGgeneral(SVGline const &eLine) {
    svgtype = SvgType::line;
    line = eLine;
  }
  SVGgeneral(SVGellipse const &eEll) {
    svgtype = SvgType::ellipse;
    ellipse = eEll;
  }
  SVGgeneral(SVGtext const &eText) {
    svgtype = SvgType::text;
    text = eText;
  }
  SvgType svgtype;
  SVGclippath clip;
  SVGpolyline polyline;
  SVGbezier bezier;
  SVGline line;
  SVGellipse ellipse;
  SVGtext text;
};

struct SVGplotDescription {
  std::vector<SVGgeneral> ListGeneral;
  /*  std::vector<SVGclippath> ListClip;
  std::vector<SVGpolyline> ListPolyline;
  std::vector<SVGbezier> ListBezier;
  std::vector<SVGline> ListLine;
  std::vector<SVGellipse> ListEllipse;*/
  double height;
  double width;
  double scale_factor;
  double add_offsetX;
  double add_offsetY;
  int FrameOption;
  int RoundMethod;
};

void GeneralWriteSVGfile(std::string const &eFile,
                         SVGplotDescription const &eSVGplot) {
  //
  // First compute the min/max of the figure.
  //
  double MinX = 0, MaxX = 0, MinY = 0, MaxY = 0;
  bool IsFirst = true;
  auto UpdateMinMaxXY = [&](coor const &pt) -> void {
    if (IsFirst) {
      MaxX = pt.x;
      MaxY = pt.y;
      MinX = pt.x;
      MinY = pt.y;
      IsFirst = false;
    } else {
      if (pt.x > MaxX)
        MaxX = pt.x;
      if (pt.x < MinX)
        MinX = pt.x;
      if (pt.y > MaxY)
        MaxY = pt.y;
      if (pt.y < MinY)
        MinY = pt.y;
    }
  };
  int nbClip = 0;
  int nbLine = 0;
  int nbPolyline = 0;
  int nbBezier = 0;
  int nbEllipse = 0;
  int nbText = 0;
  for (auto &eGeneral : eSVGplot.ListGeneral) {
    if (eGeneral.svgtype == SvgType::clip) {
      nbClip++;
    }
    if (eGeneral.svgtype == SvgType::line) {
      UpdateMinMaxXY(eGeneral.line.ePt);
      UpdateMinMaxXY(eGeneral.line.fPt);
      nbLine++;
    }
    if (eGeneral.svgtype == SvgType::polyline) {
      for (auto &eCoor : eGeneral.polyline.ListCoor)
        UpdateMinMaxXY(eCoor);
      nbPolyline++;
    }
    if (eGeneral.svgtype == SvgType::bezier) {
      UpdateMinMaxXY(eGeneral.bezier.pointM);
      UpdateMinMaxXY(eGeneral.bezier.point2);
      nbBezier++;
    }
    if (eGeneral.svgtype == SvgType::ellipse) {
      UpdateMinMaxXY(eGeneral.ellipse.c);
      nbEllipse++;
    }
    if (eGeneral.svgtype == SvgType::text) {
      UpdateMinMaxXY(eGeneral.text.point);
      nbText++;
    }
  }
  std::cerr << "nbClip = " << nbClip << "\n";
  std::cerr << "nbLine = " << nbLine << "\n";
  std::cerr << "nbPolyline = " << nbPolyline << "\n";
  std::cerr << "nbBezier = " << nbBezier << "\n";
  std::cerr << "nbEllipse = " << nbEllipse << "\n";
  std::cerr << "nbText = " << nbText << "\n";
  //
  std::cerr << "SVG: X(min/max)=" << MinX << " / " << MaxX << "\n";
  std::cerr << "SVG: Y(min/max)=" << MinY << " / " << MaxY << "\n";
  //
  // Second compute the FrameOption
  //
  double scale_factor, add_offsetX, add_offsetY;
  double height = 0, width = 0;
  bool FrameInit = false;
  // We choose the scaling factor directly from the SVGplot
  if (eSVGplot.FrameOption == 0) {
    height = eSVGplot.height;
    width = eSVGplot.width;
    scale_factor = eSVGplot.scale_factor;
    add_offsetX = eSVGplot.add_offsetX;
    add_offsetY = eSVGplot.add_offsetY;
    FrameInit = true;
  }
  // We have a fixed frame in X and Y and we adjust accordingly.
  if (eSVGplot.FrameOption == 1) {
    height = eSVGplot.height;
    width = eSVGplot.width;
    std::cerr << "height=" << height << " width=" << width << "\n";
    double FrameX = eSVGplot.width;
    double FrameY = eSVGplot.height;
    double scale_factorX = FrameY / (MaxX - MinX);
    double scale_factorY = FrameX / (MaxY - MinY);
    double MidX = (MaxX + MinX) / 2;
    double MidY = (MaxY + MinY) / 2;
    std::cerr << "MidX=" << MidX << " MidY=" << MidY << "\n";
    scale_factor = std::min(scale_factorX, scale_factorY);
    // add_offsetX + scale_factor*MidX = FrameX/2;
    add_offsetX = FrameX / 2 - scale_factor * MidX;
    add_offsetY = FrameY / 2 - scale_factor * MidY;
    FrameInit = true;
  }
  if (!FrameInit) {
    std::cerr << "FrameOption has not been used\n";
    std::cerr << "Please correct this\n";
    throw TerminalException{1};
  }
  auto GetStringValue = [&](double const &eVal,
                            double const &add_offset) -> std::string {
    double eValM = add_offset + eVal * scale_factor;
    if (eSVGplot.RoundMethod == 1)
      return DoubleTo4dot2f(eValM);
    if (eSVGplot.RoundMethod == 2)
      return DoubleToString(eValM);
    if (eSVGplot.RoundMethod == 3)
      return DoubleToString(eValM);
    std::cerr << "Failed to find relevant function\n";
    throw TerminalException{1};
  };
  auto GetStringValueX = [&](double const &eVal) -> std::string {
    return GetStringValue(eVal, add_offsetX);
  };
  auto GetStringValueY = [&](double const &eVal) -> std::string {
    return GetStringValue(eVal, add_offsetY);
  };
  auto GetStringPair = [&](coor const &pt) -> std::string {
    return GetStringValueX(pt.x) + " " + GetStringValueY(pt.y);
  };
  //
  // The attribute functionalities
  //
  auto StringColor = [&](std::vector<int> const &eV) -> std::string {
    return "rgb(" + IntToString(eV[0]) + "," + IntToString(eV[1]) + "," +
           IntToString(eV[2]) + ")";
  };
  auto f_marker = [&](std::string const &str) -> std::string {
    if (str == "")
      return "";
    return " marker-end=\"url(#" + str + ")\"";
  };
  auto f_clip = [&](std::string const &str) -> std::string {
    if (str == "")
      return "";
    return " clip-path=\"url(#" + str + ")\"";
  };
  auto GetQualityString = [&](SVGqualInfo const &eQual) -> std::string {
    return "style=\"stroke:" + StringColor(eQual.color) +
           ";stroke-width:" + DoubleToString(eQual.Size) + "\"" +
           f_marker(eQual.MarkerEnd) + f_clip(eQual.clip);
  };
  auto GetQualityStringPolyline =
      [&](SVGqualInfoPolyline const &eQual) -> std::string {
    return "style=\"fill:" + StringColor(eQual.colorfill) +
           ";stroke:" + StringColor(eQual.colorstroke) +
           ";stroke-width:" + DoubleToString(eQual.Size) + "\"" +
           f_marker(eQual.MarkerEnd) + f_clip(eQual.clip);
  };
  auto GetQualityStringEllipse = [&](SVGqualInfo const &eQual) -> std::string {
    return " style=\"fill:" + StringColor(eQual.color) + "\"" +
           f_clip(eQual.clip);
  };
  //
  // First the preamble
  //
  std::ofstream os(eFile);
  os << std::fixed;
  os << std::setprecision(9);
  os << "<svg height=\"" << height << "\" width=\"" << width << "\">\n";
  //
  // Now the major loop for the data printing
  //
  for (auto &eGeneral : eSVGplot.ListGeneral) {
    if (eGeneral.svgtype == SvgType::clip) {
      os << "  <defs>\n";
      os << "    <clipPath id=\"" << eGeneral.clip.name << "\">\n";
      os << "      <rect x=\"" << eGeneral.clip.x << "\" y=\""
         << eGeneral.clip.y << "\" width=\"" << eGeneral.clip.width
         << "\" height=\"" << eGeneral.clip.height << "\" />\n";
      os << "    </clipPath>\n";
      os << "  </defs>\n";
    }
    if (eGeneral.svgtype == SvgType::line) {
      coor ePt = eGeneral.line.ePt;
      coor fPt = eGeneral.line.fPt;
      os << "  <line x1=\"" << GetStringValueX(ePt.x) << "\" y1=\""
         << GetStringValueY(ePt.y) << "\" x2=\"" << GetStringValueX(fPt.x)
         << "\" y2=\"" << GetStringValueY(fPt.y) << "\" "
         << GetQualityString(eGeneral.line.eQual) << " />\n";
    }
    if (eGeneral.svgtype == SvgType::polyline) {
      os << "<polyline points=\"";
      bool IsFirst = true;
      for (auto &ePt : eGeneral.polyline.ListCoor) {
        if (!IsFirst)
          os << " ";
        IsFirst = false;
        os << GetStringValueX(ePt.x) << "," << GetStringValueY(ePt.y);
      }
      os << "\" " << GetQualityStringPolyline(eGeneral.polyline.eQual)
         << " />\n";
    }
    if (eGeneral.svgtype == SvgType::bezier) {
      auto eBez = eGeneral.bezier;
      os << "  <path d=\"M" << GetStringPair(eBez.pointM) << " C "
         << GetStringPair(eBez.pointC) << ", " << GetStringPair(eBez.point1)
         << ", " << GetStringPair(eBez.point2) << "\" fill=\"none\" "
         << GetQualityString(eBez.eQual) << " />\n";
    }
    if (eGeneral.svgtype == SvgType::ellipse) {
      auto eEll = eGeneral.ellipse;
      os << "  <ellipse";
      os << " cx=\"" << GetStringValueX(eEll.c.x) << "\"";
      os << " cy=\"" << GetStringValueY(eEll.c.y) << "\"";
      os << " rx=\"" << GetStringValue(eEll.r.x, 0) << "\"";
      os << " ry=\"" << GetStringValue(eEll.r.y, 0) << "\"";
      os << " " << GetQualityStringEllipse(eEll.eQual);
      os << " />\n";
    }
    if (eGeneral.svgtype == SvgType::text) {
      auto eText = eGeneral.text;
      os << "  <text";
      os << " x=\"" << GetStringValueX(eText.point.x) << "\"";
      os << " y=\"" << GetStringValueY(eText.point.y) << "\"";
      os << " font-size=\"" << eText.Size << "\"";
      os << " fill=\"" << StringColor(eText.color) << "\">";
      os << eText.string << "</text>\n";
    }
  }
  os << "</svg>\n";
}

#endif // SRC_SVG_SVGFUNCTIONS_H_
