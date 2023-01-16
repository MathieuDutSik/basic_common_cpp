// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#include "Namelist.h"

int main() {
  std::map<std::string, std::string> ListBoolValues1_doc;
  ListBoolValues1_doc["AdvancedTerminationCriterion"] = "Default: F\n\
This is about whether to used the advanced Balinski termination criterion";
  SingleBlock BlockDATA;
  BlockDATA.setListBoolValues(ListBoolValues1_doc);
}
