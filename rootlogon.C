{
   gInterpreter->AddIncludePath("include");
   gInterpreter->AddIncludePath("CppTools/include");
   gInterpreter->AddIncludePath("ROOTTools/include");
   gInterpreter->AddIncludePath("ProgressBar/include");
   gInterpreter->AddIncludePath("yaml-cpp/include");
   gInterpreter->AddIncludePath("/usr/include/jsoncpp"); // for debian
   
   gSystem->Load("CppTools/lib/libErrorHandler.so");
   gSystem->Load("CppTools/lib/libStrTools.so");
   gSystem->Load("CppTools/lib/libIOTools.so");
   gSystem->Load("CppTools/lib/libBox.so");

   gSystem->Load("ProgressBar/lib/libPBar.so");

   gSystem->Load("libjsoncpp.so");

   gSystem->Load("lib/libInputReader.so");
}
