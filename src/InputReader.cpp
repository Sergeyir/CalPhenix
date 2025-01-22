// $SOURCE$
//------------------------------------------------------------------------------------------------
//                          InputReader classes realisation
//------------------------------------------------------------------------------------------------
// InputReader
//
// ** Code for use in PHENIX related projects **
//
// Author: Sergei Antsupov
// Email: antsupov0124@gmail.com
//
/**
 * Basic classes wrapper for reading .json and .yaml files
 * This code was created to avoid for all same checks for the input file to
 * be copied in each source file in which .json and/or .yaml files are read
 **/
//------------------------------------------------------------------------------------------------

#ifndef INPUT_READER_CPP
#define INPUT_READER_CPP

#include "../include/InputReader.hpp"

InputJSONReader::InputJSONReader() {};

InputJSONReader::InputJSONReader(const std::string& inputFileOrDir, const std::string& inputType)
{
   OpenFile(inputFileOrDir, inputType);
}

void InputJSONReader::OpenFile(const std::string& inputFileOrDir, const std::string& inputType)
{
   const std::filesystem::path inputFilePath(inputFileOrDir);
   if (std::filesystem::is_directory(inputFilePath))
   {
      if (inputType == "") PrintError("InputJSONReader::OpenFile: File type was not specified");
      inputFileName = inputFileOrDir + "/" + inputType + ".json";
   }
   else
   {
      std::string fileExt = inputFileOrDir;
      fileExt.erase(0, inputFileOrDir.size() - 5);
      
      if (fileExt != ".json") 
      {
         PrintError("InputJSONReader::OpenFile: Input file must have .json extention.");
      }
      inputFileName = inputFileOrDir;
   }

   CheckInputFile(inputFileName);
   
   inputFile.open(inputFileName.c_str(), std::ifstream::binary);
   inputFile >> inputFileContents;
}

void InputJSONReader::CheckStatus(const std::string& status)
{
   if (inputFileContents["status"].asString() != status)
   {
      PrintError("Status mismatch in file " + inputFileName + ": " +
                 inputFileContents["status"].asString() + " vs " + status);
   }
}

Json::Value InputJSONReader::operator[](const std::string& field)
{
   return inputFileContents[field];
}

InputJSONReader::~InputJSONReader() {};

#endif /* INPUT_READER_CPP */
