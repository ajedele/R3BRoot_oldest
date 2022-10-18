/******************************************************************************
 *   Copyright (C) 2019 GSI Helmholtzzentrum für Schwerionenforschung GmbH    *
 *   Copyright (C) 2019 Members of R3B Collaboration                          *
 *                                                                            *
 *             This software is distributed under the terms of the            *
 *                 GNU General Public Licence (GPL) version 3,                *
 *                    copied verbatim in the file "LICENSE".                  *
 *                                                                            *
 * In applying this license GSI does not waive the privileges and immunities  *
 * granted to it by virtue of its status as an Intergovernmental Organization *
 * or submit itself to any jurisdiction.                                      *
 ******************************************************************************/

//*-- AUTHOR : D. Kresan
//*-- Created : 18/05/2015

/////////////////////////////////////////////////////////////
//
//  R3BPspxContFact
//
//  Factory for the parameter containers in libR3BPspx
//
/////////////////////////////////////////////////////////////

#include "R3BPspxContFact.h"

#include "R3BGeoPspxPar.h"
#include "R3BTGeoPar.h"

#include "FairParAsciiFileIo.h"
#include "FairParRootFileIo.h"
#include "FairRuntimeDb.h"

#include "R3BPspxCalPar.h"
#include "R3BPspxHitPar.h"

#include "TClass.h"

#include <iomanip>
#include <iostream>


static R3BPspxContFact gR3BPspxContFact;

R3BPspxContFact::R3BPspxContFact()
{
    // Constructor (called when the library is loaded)
    fName = "R3BPspxContFact";
    fTitle = "Factory for parameter containers in libR3BPspx";
    setAllContainers();
    FairRuntimeDb::instance()->addContFactory(this);
}

void R3BPspxContFact::setAllContainers()
{
     FairContainer* p1= new FairContainer("PspxCalPar", "Pspx Cal Parameters", "PspxCalParContext");
     p1->addContext("PspxCalParContext");
     containers->Add(p1);

     FairContainer* p2= new FairContainer("PspxHitPar", "Pspx Hit Parameters", "PspxHitParContext");
     p2->addContext("PspxHitParContext");
     containers->Add(p2);


     FairContainer* p3 = new FairContainer("PspxGeoPar", "Pspx geometry parameters", "GeometryParameterContext");
     p3->addContext("GeometryParameterContext");
     containers->Add(p3);
}

FairParSet* R3BPspxContFact::createContainer(FairContainer* c)
{
    /** Pspxls the constructor of the corresponding parameter container.
     * For an actual context, which is not an empty string and not the default context
     * of this container, the name is concatinated with the context. */

    const char* name = c->GetName();
    cout << " R3BPspxContFact: Create container name " << name << endl;
    FairParSet* p = 0;
    if (strcmp(name,"PspxCalPar")==0) 
    {
      p=new R3BPspxCalPar(c->getConcatName().Data(),c->GetTitle(),c->getContext());
      cout << " Inside CalPar Container creator " << name << endl;
    }
    else if (strcmp(name,"PspxHitPar")==0) 
    {
      p=new R3BPspxHitPar(c->getConcatName().Data(),c->GetTitle(),c->getContext());
    }

    else if (strcmp(name, "PspxGeoPar") == 0)
    {
        p = new R3BTGeoPar(c->getConcatName().Data(), c->GetTitle(), c->getContext());
    }

    return p;
}

void R3BPspxContFact::activateParIo(FairParIo* io)
{
	// activated the input/output class for the parameters
	// needed by the Sts
}

ClassImp(R3BPspxContFact);
