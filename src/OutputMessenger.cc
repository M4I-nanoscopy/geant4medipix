#include "RunAction.hh"
#include "G4UImessenger.hh"
#include "G4UIcmdWithAString.hh"
#include "OutputMessenger.hh"
#include "G4UserRunAction.hh"
#include "G4RunManager.hh"

OutputMessenger::OutputMessenger()
:runAction((RunAction *) G4RunManager::GetRunManager()->GetUserRunAction())
{
	dir = new G4UIdirectory("/run/");
	dir->SetGuidance("Change the output file name.");

	nameCmd = new G4UIcmdWithAString("/run/output",this);
	nameCmd->SetGuidance("Change the output file name.");
	nameCmd->SetGuidance(" (output.hdf5 is default)");
	nameCmd->SetParameterName("fileName",true,true);
	nameCmd->SetDefaultValue("output.hdf5");
	nameCmd->AvailableForStates(G4State_PreInit,G4State_Idle,G4State_Init);
}

OutputMessenger::~OutputMessenger()
{
  delete nameCmd;
  delete dir;
}

void OutputMessenger::SetNewValue(
  G4UIcommand * command,G4String newValues)
{
	if( command==nameCmd )
	{
		if(newValues.strip() != "")
		{
			runAction->SetName(newValues);
		}
	}
}

G4String OutputMessenger::GetCurrentValue(G4UIcommand * command)
{
  G4String cv;
  if( command==nameCmd )
  {
	  cv = runAction->GetName();
  }
  return cv;
}