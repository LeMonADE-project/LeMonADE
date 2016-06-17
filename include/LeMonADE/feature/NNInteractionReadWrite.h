#ifndef NNINTERACTIONREADWRITE_H
#define NNINTERACTIONREADWRITE_H

#include<iostream>

/******************************************************************************
 *definition of read and write classes used by FeatureNNInteraction
 *and FeatureNNInteractionPartialEV
 ******************************************************************************/
template < class IngredientsType>
class ReadInteraction: public ReadToDestination<IngredientsType>
{
public:
    ReadInteraction(IngredientsType& i):ReadToDestination<IngredientsType>(i){}
    virtual ~ReadInteraction(){}
    virtual void execute();
};


template <class IngredientsType>
class WriteInteraction:public AbstractWrite<IngredientsType>
{
public:
    WriteInteraction(const IngredientsType& i)
        :AbstractWrite<IngredientsType>(i){this->setHeaderOnly(true);}

    virtual ~WriteInteraction(){}

    virtual void writeStream(std::ostream& strm);
};

///////////////////////////////////////////////////////////////////////////////
template<class IngredientsType>
void ReadInteraction<IngredientsType>::execute()
{
    IngredientsType& ingredients=this->getDestination();
    std::istream& file=this->getInputStream();


    uint32_t typeA,typeB;
    double interactionConstant;


    file>>typeA;

    if(file.fail())
    {
        std::stringstream messagestream;
        messagestream<<"ReadInteraction::execute()\n"
                    <<"Could not read first type\n";
        throw std::runtime_error(messagestream.str());
    }

    //now read second type
    file>>typeB;

    if(file.fail())
    {
        std::stringstream messagestream;
        messagestream<<"ReadInteraction::execute()\n"
                    <<"Could not read second type\n";
        throw std::runtime_error(messagestream.str());
    }

    //now read interactionConstant
    file>>interactionConstant;

    if(file.fail())
    {
        std::stringstream messagestream;
        messagestream<<"ReadInteraction::execute()\n"
                    <<"Could not read interaction constant.\n";
        throw std::runtime_error(messagestream.str());
    }


    //now save the interaction tuple just read from the file
    ingredients.setInteraction(typeA,typeB,interactionConstant);


}



template<class IngredientsType>
void WriteInteraction<IngredientsType>::writeStream(std::ostream& stream)
{
    size_t nSpecies=10;
    stream<<"## nearest neighbor interactions between types in kT (default 0.0kT)\n";

    for(size_t typeA=1;typeA<=nSpecies;typeA++)
    {
        for(size_t typeB=1;typeB<typeA;typeB++)
        {
            if(this->getSource().getInteraction(typeA,typeB)!=0.0)
            {
                stream<<"#!nn_interaction "<<typeB<<" "<<typeA<<" "<<this->getSource().getInteraction(typeB,typeA)<<"\n";
            }

        }

    }
    stream<<"\n\n";

}




#endif // NNINTERACTIONREADWRITE_H
