/** 
 *  @file   InputYAMLReader.hpp 
 *  @brief  Contains declaration of class InputYAMLReader
 *
 *  This file is a part of a project CalPhenix (https://github.com/Sergeyir/CalPhenix).
 *
 *  @author Sergei Antsupov (antsupov0124@gmail.com)
 **/
#ifndef INPUT_READER_HPP
#define INPUT_READER_HPP

#include <string>
#include <filesystem>

#include "yaml-cpp/yaml.h"

#include "ErrorHandler.hpp"

/*! @class InputYAMLReader
 * @brief Class InputYAMLReader can be used to simplify the work with .yaml input files with the same formatting
 *
 * The same formatting includes.yaml field "status" which shortly describes what the .yaml file is for. With this class you don't need to check the existence of files since InputYAMLReader does it automaticaly. Also this class is useful since in the code you can just point to InputYAMLReader object what file it is expected to read and then pass only directory to it (see InputYAMLReader constructor with parameters). This simplifies the process and saves time typing the name of the file if you pass the directory instead of directory + input file name into compiled executable call that expects input file name or directory as an argument.
 */
class InputYAMLReader
{
   public:
   
   ///@brief Default constructor
   InputYAMLReader();
   /*! @brief Constructor with parameters
    * See InputYAMLReader::Open(const std::string& inputFileOrDir, const std::string& inputType = "") for details on parameters and checks
    */
   InputYAMLReader(const std::string& inputFileOrDir, const std::string& inputType = "");
   /*! @brief Opens the .yaml file
    * @param[in] inputFileOrDir input file name or the directory the input file is in. If inputFileOrDir is not directory InputYAMLReader will try to open this file unless it is not .yaml file in which case InputYAMLReader will print error and will exit the program with exit code 1.
    * @param[in] inputType input type of the file (e.g..yaml field "status" in the file). If inputFileOrDir parameter is a directory name InputYAMLReader will try to open file inputFileOrDir + "/" + inputType.yaml.
    * If the input file does not contain field "status" with value equal to inputType InputYAMLReader will print error and exit the program with exit code 1.
    */
   void OpenFile(const std::string& inputFileOrDir, const std::string& inputType = "");
   /*! @brief Perfoms the check for .yaml field "status"
    * If the input file does not contain field "status" with value equal to passed status value InputYAMLReader will print error and exit the program with exit code 1.
    */
   void CheckStatus(const std::string& status);
   /// @brief Public access to YAML::Node operator[] associated with the file contents opened with InputYAMLReader object
   YAML::Node operator[](const std::string& field);
   /// @brief Default destructor
   virtual ~InputYAMLReader();

   private:
   
   /// Name of the .yaml file
   std::string inputFileName;
   /// Input file contents
   YAML::Node inputFileContents;
};

#endif /* INPUT_READER_HPP */
