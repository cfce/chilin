import os
import copy
from copy import deepcopy
import subprocess
from time import strftime, localtime
from samflow.helper import print_command_details

class AbstractCommand(object):
    def __init__(self, template=None, tool=None, param={}, input=[], output=[], name=""):
        self.name = name
        self.input = input
        self.output = output
        self.param = param
        self.template = template
        self.tool = tool
        self.result = None
        self.verbose_level = 1
        self.dry_run_mode = False
        self.resume = False
        self.fetch_output = False

        self.allow_fail = False
        self.allow_dangling = False

        # if the parent of a command is itself, it's a root
        self._parent = self
        self._commands = []
        self._allow_zero_byte_file = True

    def add_back(self, command):
        """ For workflow only, add a command into current workflow """
        raise NotImplemented

    def add_front(self, command):
        raise NotImplemented

    def update(self, **kwargs):
        for k, v in kwargs.items():
            getattr(self, k).update(v)
        return self

    def set(self, **kwargs):
        for k, v in kwargs.items():
            setattr(self, k, v)
        return self

    @property
    def have_dangling(self):
        dangling_inputs = self._dangling_inputs

        if dangling_inputs:
            print(self.allow_dangling)
            if self.allow_dangling:
                self._print_log("Warn", "Dangling inputs ", dangling_inputs)
            else:
                self._print_log("Error!", "Dangling inputs ", dangling_inputs)
                return True
        return False

    def invoke(self):
        """
        Invoke the command, return True if no trouble encountered.
        Dry run mode check files `dangling` for only input
        Non-dry run mode check `missing` for both input and output
        """
        try:
            if self.have_dangling:
                return False

            if self.dry_run_mode:
                self.result = self._simulate()
                return True

            missing_i = self._missing_inputs

            if missing_i:
                self._print_log("Error!", "Missing inputs", missing_i)
                if self.allow_fail:
                    return True ## allow missing inputs
                else:
                    return False

            execute_success = self._execute()
            if self.allow_fail:
                return True
            if not execute_success:
                return False

            missing_o = self._missing_outputs
            if missing_o:
                return False
            return True
        except:
            print("Exception encountered @", self.name)
            print_command_details(self)
            raise

    def __deepcopy__(self, visit):
        return type(self)(template=deepcopy(self.template), tool=deepcopy(self.tool), param=deepcopy(self.param),
            input=deepcopy(self.input), output=deepcopy(self.output), name=deepcopy(self.name))

    def set_option(self, **args):
        """
        :type verbose_level: int
        verbose_level:  1 - only show fatal errors          (quiet mode)
                        2 - show fatal errors and warnings  (normal mode)
                        3 - show workflow details           (verbose mode)
                        4 - debug mode                      (debug mode)
        """
        for k, v in args.items():
            setattr(self, k, v)

    @property
    def clone(self):
        return copy.deepcopy(self)

    @property
    def _is_root(self):
        return self._parent == self

    def _simulate(self):
        """ Hook method for `invoke` in dry run mode: Pretend to run but not invoke anything """
        pass

    def _print_log(self, head, *args):
        if isinstance(self, ShellCommand) and head in ["Run", "Dry-run"]:
            start_with = ""
        else:
            start_with = "#"
        print start_with, head, self.name, "[", strftime("%Y-%m-%d %H:%M:%S", localtime()), "]: \n", " ".join(map(str, args))

    def _execute(self):
        """
        Hook method for `invoke`
        Return True if no trouble encountered
        """
        return True

    @property
    def _dummy_files(self):
        """
        Return a list of files that produced by leaves before current node in the tree
        """
        if self._is_root:
            # if current command is not a leaf of a tree, it shouldn't have dummy files.
            return []
        ret = []
        for a_command in self._root:
            if a_command == self:
                break
            else:
                ret += a_command._outputs
        return ret

    def _missing(self, files):
        """ use step name to identify error step,
        self.name """
        missing = []

        try:
            for i in files:
                if not os.path.exists(i):
                    missing.append(i)
                    continue
                if self._allow_zero_byte_file: ## default: accept 0 byte file
                    continue
                if os.path.isfile(i) and os.path.getsize(i) == 0:
                    missing.append(i)
                    continue
            return missing
        except:
            print("Exception encountered @", self.name, self.template)
            raise

    @property
    def _missing_inputs(self):
        """
        Return a list of files that are current command's input but doesn't exist in filesystem .
        This method is called before real invoke current command.
        """
        return self._missing(self._inputs)

    @property
    def _missing_outputs(self):
        """
        Return a list of files that are current command's output but doesn't exist in filesystem .
        This method is called after real invoke current command.
        """
        return self._missing(self._outputs)

    @property
    def _dangling_inputs(self):
        """
        Return a list of files that are current commands' input but:
        (1) doesn't exist in filesystem
        (2) can't be found as some commands' output before current command

        This Hook method is called:
        (1) on each leaf before both dry run and real run

        If current command doesn't belong to a tree, just return missing inputs
        """
        if self._is_root:
            return self._missing_inputs
        else:
            return [i for i in self._missing_inputs if i not in self._dummy_files]

    @property
    def _root(self):
        """
        Return the root of the tree in which current command is
        """
        return self if self._is_root else self._parent._root

    def _collect(self, obj):
        ret = []
        if isinstance(obj, dict):
            for i in obj.values():
                ret.extend(self._collect(i))
        elif isinstance(obj, str):
            ret = [obj]
        elif isinstance(obj, list):
            for i in obj:
                ret.extend(self._collect(i))

        return ret

    @property
    def _inputs(self):
        """ Return the inputs as a list """
        return self._collect(self.input)

    @property
    def _outputs(self):
        """ Return the outputs as a list """
        return self._collect(self.output)


class ShellCommand(AbstractCommand):
    def __init__(self, template=None, tool=None, param={}, input=[], output=[], name=""):
        AbstractCommand.__init__(self, template, tool, param, input, output, name)

        self.fetch_output = False

    def _simulate(self):
        self._print_log("Dry-run", self._render())

    def _execute(self):
        cmd_rendered = self._render()
        can_skip = not self._missing_outputs
        if self.resume:
            if can_skip:
                #print("Resumed from existing result.. Skip")
                return True
            else:
                self._print_log("Run", cmd_rendered)
        else:
            self._print_log("Run", cmd_rendered)
        if self.fetch_output:
            try:
                self.result = subprocess.check_output(cmd_rendered, shell=True, universal_newlines=True,
                    executable="/bin/bash")
            except subprocess.CalledProcessError:
                return False
        else:
            try:
                self.result = subprocess.check_call(cmd_rendered, shell=True, executable="/bin/bash")
            except subprocess.CalledProcessError:
                return False
        return True

    def _render(self):
        """ Method that return the rendered content  """
        cmd = self.template.format(input=self.input, output=self.output, param=self.param, tool=self.tool)
        return cmd

    def set_stdout_collecting(self):
        self.fetch_output = True
        return self

    @property
    def _have_dangling_tools(self):
        if not self.tool:
            return False
        elif which(self.tool):
            return False
        else:
            return True

    @property
    def have_dangling(self):
        have_dangling_inputs = AbstractCommand.have_dangling
        have_dangling_tool = self._have_dangling_tools
        have_dangling = have_dangling_inputs and have_dangling_tool

        if have_dangling_tool:
            # print self.allow_dangling
            if self.allow_dangling:
                self._print_log("Warn", "Dangling tool ", self.tool)
                return have_dangling
            else:
                self._print_log("Error!", "Dangling tool", self.tool)
                return True
        return have_dangling


class PythonCommand(AbstractCommand):
    def __init__(self, template=None, tool=None, param={}, input=[], output=[], name=""):
        AbstractCommand.__init__(self, template, tool, param, input, output, name)

    def _render(self):
        return "%s < %s > %s" % (self.template, self._inputs, self._outputs)

    def _simulate(self):
        self._print_log("Dry-run", self._render())
        return None

    def _execute(self):
        cmd_rendered = self._render()
        can_skip = not self._missing_outputs ## output test
        if self.resume:
            if can_skip:
                #print("Resumed from existing result.. Skip") ## suppress information
                return True
            else:
                self._print_log("Run", cmd_rendered)
        else:
            self._print_log("Run", cmd_rendered)
        try:
            self.template(input=self.input, output=self.output, param=self.param) ## in case python internal error
        except Exception as e:
            print e
        return True


def which(program):
    import os

    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = os.path.expanduser(path.strip('"'))
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return None
