// $HEADER$
//------------------------------------------------------------------------------------------------
//                          InputReader classes declaration
//------------------------------------------------------------------------------------------------
// InputReader
//
// ** Code for use in PHENIX related projects **
//
// Author: Sergei Antsupov
// Email: antsupov0124@gmail.com
//
/**
 * Basic classes wrappers for reading .json and .yaml files
 * This code was created to avoid all checks for the input fiel to
 * be copied in each source file where .json and/or .yaml files are read
 **/
//------------------------------------------------------------------------------------------------

#ifndef INPUT_READER_HPP
#define INPUT_READER_HPP

#include <string>
#include <fstream>
#include <filesystem>

#include "json/json.h"
#include "json/value.h"

#include "ErrorHandler.hpp"

class InputJSONReader
{
   public:
   
   InputJSONReader();
   InputJSONReader(const std::string& inputFileOrDir, const std::string& inputType = "");
   void OpenFile(const std::string& inputFileOrDir, const std::string& inputType = "");
   void CheckStatus(const std::string& status);
   Json::Value operator[](const std::string& field);
   ~InputJSONReader();

   private:
   
   std::string inputFileName;
   std::ifstream inputFile;
   Json::Value inputFileContents;
};

#endif /* INPUT_READER_HPP */
