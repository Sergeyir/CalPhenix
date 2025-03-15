/** 
 *  @file   InputYAMLReader.cpp 
 *  @brief  Contains declaration of class InputYAMLReader
 *
 *  In order to use InputYAMLReader class libraries libInputReader.so, libErrorHandler.so must be loaded
 *
 *  This file is a part of a project CalPhenix (https://github.com/Sergeyir/CalPhenix).
 *
 *  @author Sergei Antsupov (antsupov0124@gmail.com)
 **/
#ifndef INPUT_READER_CPP
#define INPUT_READER_CPP

#include "../include/InputYAMLReader.hpp"

InputYAMLReader::InputYAMLReader() {};

InputYAMLReader::InputYAMLReader(const std::string& inputFileOrDir, const std::string& inputType)
{
   OpenFile(inputFileOrDir, inputType);
}

void InputYAMLReader::OpenFile(const std::string& inputFileOrDir, const std::string& inputType)
{
   const std::filesystem::path inputFilePath(inputFileOrDir);
   if (std::filesystem::is_directory(inputFilePath))
   {
      if (inputType == "") CppTools::PrintError("InputYAMLReader::OpenFile: File type was not specified");
      inputFileName = inputFileOrDir + "/" + inputType + ".yaml";
   }
   else
   {
      std::string fileExt = inputFileOrDir;
      fileExt.erase(0, inputFileOrDir.size() - 5);
      
      if (fileExt != ".yaml") 
      {
         CppTools::PrintError("InputYAMLReader::OpenFile: Input file must have .yaml extention while " + fileExt + " type file was provided");
      }
      inputFileName = inputFileOrDir;
   }

   CppTools::CheckInputFile(inputFileName);
   
   inputFileContents = YAML::LoadFile(inputFileName);
}

void InputYAMLReader::CheckStatus(const std::string& status)
{
   if (inputFileContents["status"].as<std::string>() != status)
   {
      CppTools::PrintError("Status mismatch in file " + inputFileName + ": " +
                           inputFileContents["status"].as<std::string>() + " vs " + status);
   }
}

YAML::Node InputYAMLReader::operator[](const std::string& field)
{
   return inputFileContents[field];
}

InputYAMLReader::~InputYAMLReader() {};

#endif /* INPUT_READER_CPP */
