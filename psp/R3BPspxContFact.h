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

#ifndef R3BPSPXCONTFACT_H
#define R3BPSPXCONTFACT_H

#include "FairContFact.h"

class FairContainer;

class R3BPspxContFact : public FairContFact
{
  private:
    void setAllContainers();

  public:
    R3BPspxContFact();
    ~R3BPspxContFact() {}
    FairParSet* createContainer(FairContainer*);
    void activateParIo(FairParIo* io);
    ClassDef(R3BPspxContFact, 0) // Factory for all PSPX parameter containers
};

#endif /* !R3BPSPXCONTFACT_H */
