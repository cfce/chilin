from samflow.command import AbstractCommand
from samflow.helper import print_command_details


class Workflow(AbstractCommand):
    def __init__(self, template=None, tool=None, param = {},input=[], output=[], name=""):
        AbstractCommand.__init__(self, template=None, tool=None, param = {},input=[], output=[], name=name)
        self._commands = []

    def __iter__(self):
        for cmd_i in self._commands:
            if isinstance(cmd_i, Workflow):
                for cmd_i_j in cmd_i:
                    yield cmd_i_j
            else:
                yield cmd_i

    def add_back(self, command):
        command._parent = self
        self._commands.append(command)
        return self

    def add_front(self, command):
        command._parent = self
        self._commands.insert(0, command)
        return self

    @property
    def _dangling_inputs(self):
        dangling_dict = {}
        for cmd in self:
            not_dangling = not cmd._dangling_inputs
            if not_dangling:
                continue
            if cmd.name in dangling_dict:
                # update by union of two set
                dangling_dict[cmd.name] |= set(cmd._dangling_inputs)
            else:
                dangling_dict[cmd.name] = set(cmd._dangling_inputs)
        return dangling_dict

    def _have_render_error(self):
        error_keys = []
        for cmd in self:
            try:
                cmd._render()
            except KeyError as ke:
                error_keys.append([cmd.name, ke])
            except:
                print("Exception encountered when rendering @", cmd.name)
                print_command_details(cmd)
                raise

        if error_keys:
            self._print_log("KeysError!", error_keys)
            return True
        return False

    def invoke(self):
        if self._have_render_error(): ## template rendering test
            return False

        if self.have_dangling:  ## dangling input test for dry run check
            return False

        for cmd in self:
            success_invoked = cmd.invoke()
            if not success_invoked:
                print(cmd.name, cmd._render())
                print("{0:!^80}".format("Error happened! Workflow stopped!"))
                return False
        print("{0:-^80}".format("Workflow finished successfully"))
        return True

    def set_option(self, **args):
        AbstractCommand.set_option(self, **args)
        for cmd in self._commands:
            cmd.set_option(**args)

def attach_back(Workflow, AbstractCommand):
    Workflow.add_back(AbstractCommand)
    return AbstractCommand

def attach_front(Workflow, AbstractCommand):
    Workflow.add_front(AbstractCommand)
    return AbstractCommand
