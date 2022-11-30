from command_scriptor import scriptorSwitcher
from constants import Constants
from filenames import FilenameList
from base import BaseCommander
from base_extension import BaseCommanderExtended

class Optimizer:
    scriptor = None
    const = None
    fl = None
    base = None

    def __init__(self, config):
        scriptorType = config.getInput(['command_scriptor', 'type'], datatype=str)
        self.scriptor = scriptorSwitcher.get(scriptorType)

        self.const = Constants(config)
        self.fl = FilenameList(config, self.const)

        if (config.getInput(['time_splitting', 'use_state_mollifier'], fallback=False)):
            self.base = BaseCommanderExtended(config, self.scriptor, self.const, self.fl)
        else:
            self.base = BaseCommander(config, self.scriptor, self.const, self.fl)
        return
