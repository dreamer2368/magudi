from command_scriptor import scriptorSwitcher
from constants import Constants

class Optimizer:
    scriptor = None
    const = None

    def __init__(self, config):
        scriptorType = config.getInput(['command_scriptor', 'type'], datatype=str)
        self.scriptor = scriptorSwitcher.get(scriptorType)

        self.const = Constants(config)
        return
