//
// Created by lucas on 31/07/17.
//

#ifndef G4MEDIPIX_OUTPUTMESSENGER_HH
#define G4MEDIPIX_OUTPUTMESSENGER_HH

#include <G4UImessenger.hh>

class G4UIcommand;
class G4UIcmdWithAString;


class RunAction;

class OutputMessenger : public G4UImessenger
{
	public:
    	OutputMessenger();
    	~OutputMessenger();

	public:
	    void SetNewValue(G4UIcommand * command,G4String newValues);
	    G4String GetCurrentValue(G4UIcommand * command);

	private:
	    RunAction * runAction;

	private:
	    G4UIdirectory * dir;
	    G4UIcmdWithAString * nameCmd;
};


#endif //G4MEDIPIX_OUTPUTMESSENGER_HH
