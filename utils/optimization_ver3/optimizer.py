from command_scriptor import scriptorSwitcher
from constants import Constants
from filenames import FilenameList

class Optimizer:
    scriptor = None
    const = None
    fl = None

    def __init__(self, config):
        scriptorType = config.getInput(['command_scriptor', 'type'], datatype=str)
        self.scriptor = scriptorSwitcher.get(scriptorType)

        self.const = Constants(config)
        self.fl = FilenameList(config, self.const)
        return
