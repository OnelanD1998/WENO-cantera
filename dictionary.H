#include <fstream>
#include <string>
#include <sstream>

#ifndef DICTIONARY_H
#define DICTIONARY_H

class dictionary
{
private:
    std::ifstream fileDictionary;
    std::string dictionaryName_;
public:
    dictionary();
    dictionary(char fileName[]);
    dictionary(std::string fileName);

    template<typename TT>
    TT lookUp(std::string headName)//char headName[])
    {
        fileDictionary.seekg(0, std::ios::beg);
        using namespace std;
        string meshDicLine;

        while(fileDictionary >> meshDicLine)
        {

            if(fileDictionary.eof()) 
            {
                std::cerr << "|******* No " << headName << " is find in " << dictionaryName_ <<"*******| \n ";
                throw std::invalid_argument("invalid_argument");
            }

            if(meshDicLine == headName)
            {   
                string temp1;
                TT tempTT;
                fileDictionary >> temp1;

                if(temp1 == "\n") throw std::invalid_argument("invalid_argument");

                std::istringstream iss(std::move(temp1));
                if(!(iss >> tempTT))
                    throw std::invalid_argument("invalid_argument");

                return tempTT;
            }
        }
    }

    ~dictionary();
};

dictionary::dictionary()
{}

dictionary::dictionary(char fileName[])
:
fileDictionary(fileName),
dictionaryName_(fileName)
{
    if(!fileDictionary.is_open()){ 
        std::cout << "the file {" << dictionaryName_ << "} is not find" <<std::endl;
        throw std::invalid_argument("invalid_argument");}
}

dictionary::dictionary(std::string fileName)
:
fileDictionary(fileName),
dictionaryName_(fileName)
{
    if(!fileDictionary.is_open()){ 
        std::cout << "the file {" << dictionaryName_ << "} is not find" <<std::endl;
        throw std::invalid_argument("invalid_argument");}
}

dictionary::~dictionary()
{
}

#endif