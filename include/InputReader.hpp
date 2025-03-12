/** 
 *  @file   InputReader.hpp 
 *  @brief  Contains declaration of class InputJSONReader
 *
 *  In order to use InputJSONReader class libraries libInputReader.so, libErrorHandler.so must be loaded
 *
 *  This file is a part of a project CalPhenix (https://github.com/Sergeyir/CalPhenix).
 *
 *  @author Sergei Antsupov (antsupov0124@gmail.com)
 **/
#ifndef INPUT_READER_HPP
#define INPUT_READER_HPP

#include <string>
#include <fstream>
#include <filesystem>

#include "json/json.h"
#include "json/value.h"

#include "ErrorHandler.hpp"

/*! @class InputJSONReader
 * @brief Class InputJSONReader can be used to simplify the work with .json input files with the same formatting
 *
 * The same formatting includes json field "status" which shortly describes what the .json file is for. With this class you don't need to check the existence of files since InputJSONReader does it automaticaly. Also this class is useful since in the code you can just point to InputJSONReader object what file it is expected to read and then pass only directory to it (see InputJSONReader constructor with parameters). This simplifies the process and saves time typing the name of the file if you pass the directory instead of directory + input file name into compiled executable call that expects input file name or directory as an argument.
 */
class InputJSONReader
{
   public:
   
   ///@brief Default constructor
   InputJSONReader();
   /*! @brief Constructor with parameters
    * See InputJSONReader::Open(const std::string& inputFileOrDir, const std::string& inputType = "") for details on parameters and checks
    */
   InputJSONReader(const std::string& inputFileOrDir, const std::string& inputType = "");
   /*! @brief Opens the .json file
    * @param[in] inputFileOrDir input file name or the directory the input file is in. If inputFileOrDir is not directory InputJSONReader will try to open this file unless it is not .json file in which case InputJSONReader will print error and will exit the program with exit code 1.
    * @param[in] inputType input type of the file (e.g. json field "status" in the file). If inputFileOrDir parameter is a directory name InputJSONReader will try to open file inputFileOrDir + "/" + inputType.json.
    * If the input file does not contain field "status" with value equal to inputType InputJSONReader will print error and exit the program with exit code 1.
    */
   void OpenFile(const std::string& inputFileOrDir, const std::string& inputType = "");
   /*! @brief Perfoms the check for .json field "status"
    * If the input file does not contain field "status" with value equal to passed status value InputJSONReader will print error and exit the program with exit code 1.
    */
   void CheckStatus(const std::string& status);
   /// @brief Public access to Json::Value operator[] associated with the file contents opened with InputJSONReader object
   Json::Value operator[](const std::string& field);
   /// @brief Default destructor
   virtual ~InputJSONReader();

   private:
   
   /// Name of the .json file
   std::string inputFileName;
   /// Input file that is opened with the object of this class
   std::ifstream inputFile;
   /// Input file contents
   Json::Value inputFileContents;
};

#endif /* INPUT_READER_HPP */
