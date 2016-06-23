from ConfigParser import ConfigParser, NoSectionError, NoOptionError
import os
import logging

class NoTreatmentData(Exception):
    pass

class ChiLinConfig(ConfigParser):
    def __init__(self, conf, args):
        ConfigParser.__init__(self)
        self._verbose_level = 1
        self._conf = ConfigParser()
        if not os.path.exists(conf):
            raise IOError("No such config file: %s" % repr(conf))
        self._conf.read(conf)
        self.root_dir = os.path.dirname(conf)
        self.pe = args.pe
        self.long = False

    def write(self, fileobj):
        self._conf.write(fileobj)

    def set_option(self, verbose_level=1):
        """
        :type verbose_level: int
        verbose_level:  1 - only show fatal errors          (quiet mode)
                        2 - show fatal errors and warnings  (normal mode)
                        3 - show workflow details           (verbose mode)
                        4 - debug mode                      (debug mode)
        """
        self._verbose_level = verbose_level

    def get(self, section, option, default = None):
        try:
            return self._conf.get(section, option)
        except NoOptionError:
            if default:
                return default
            else:
                raise

    @property
    def treatment_pairs_pe(self):
        return list(zip(self.treatment_raws, self.treatment_pair_targets["pairs"]))

    def set(self, section, option, value):
        try:
            self._conf.set(section, option, value)
        except NoOptionError:
            raise

    def get_path(self, section, option):
        return self.to_abs_path(self.get(section, option))

    def items(self, section):
        try:
            return self._conf.items(section)
        except NoSectionError:
            if self._verbose_level >= 2:
                print("Warning: No such section: ", section)
                print("This will return a empty dict")
            return {}

    @property
    def id(self):
        return self.get("basics", "id")

    @property
    def target_dir(self):
        return self.get("basics", "output")

    @property
    def prefix(self):
        return os.path.join(self.target_dir, self.id)

    @property
    def json_prefix(self):
        return os.path.join(self.category("json"), self.id)

    @property
    def latex_prefix(self):
        return os.path.join(self.category("latex"), self.id)

    @property
    def treatment_pairs(self):
        """
        one to one in single end mode,  original, target
        two to one in pair end mode [original_pair1, original_pair2], target
        """
        if not self.pe:
            return list(zip(self.treatment_raws, self.treatment_targets))
        else:
            return self.treatment_pairs_pe

    @property
    def control_pairs(self):
        if not self.pe:
            return list(zip(self.control_raws, self.control_targets))
        else:
            return self.control_pairs_pe

    @property
    def sample_pairs(self):
        return self.treatment_pairs + self.control_pairs

    def to_abs_path(self, path):
        abs_path = path
        if not os.path.isabs(path):
            #abs_path = os.path.join(self.root_dir, abs_path)
            abs_path = os.path.abspath(abs_path)
        return abs_path

    @property
    def control_raws(self):
        if self.get("basics","cont").strip():
            if not self.pe:
                return [self.to_abs_path(i.strip()) for i in self.get("basics", "cont").split(",")]
            else:
                data_list = []
                for i in self.get("basics", "cont").split(";"):
                    data_list.append([ self.to_abs_path(j.strip()) for j in i.split(",") ])
                return data_list
        return []


    @property
    def treatment_raws(self):
        """
        single end data separate by , for replicates
        pair end data separate by ; for replicates , for pairs
        """
        if self.get("basics", "treat"):
            if not self.pe:
                return [self.to_abs_path(i.strip()) for i in self.get("basics", "treat").split(",")]
            else:
                data_list = []
                for i in self.get("basics", "treat").split(";"):
                    data_list.append([ self.to_abs_path(j.strip()) for j in i.split(",") ])
                return data_list
        else:
            raise NoTreatmentData

    # previous interface only for SE
    # @property
    # def treatment_raws(self):
    #     if self.get("basics", "treat").strip():
    #         return [self.to_abs_path(i.strip()) for i in self.get("basics", "treat").split(",")]
    #     else:
    #         raise NoTreatmentData

    # previous interface only for SE
    # @property
    # def treatment_targets(self):
    #     return [os.path.join(self.target_dir,
    #         self.id + "_treat_rep" + str(num+1)
    #     ) for num in range(len(self.treatment_raws))]

    @property
    def pe(self):
        return self._pe

    @pe.setter
    def pe(self, value):
        '''setting Pair End state, True for PE, False for SE
        '''
        self._pe = value


    @property
    def treatment_targets(self):
        if not self.pe:
            return self.treatment_single_targets
        else:
            return self.treatment_pair_targets["reps"]

    @property
    def control_pairs_pe(self):
        return list(zip(self.control_raws, self.control_pair_targets["pairs"]))

    @property
    def treatment_pair_data(self):
        return self.treatment_pair_targets["pairs"]

    @property
    def treatment_single_targets(self):
        return [os.path.join(self.target_dir,
            self.id + "_treat_rep" + str(num+1)
        ) for num in range(len(self.treatment_raws))]

    @property
    def treatment_pair_targets(self):
        '''pairs: for [[rep1_pair1, rep1_pair2]],
        usually for evaluating read quality
        reps: for [rep1, rep2],
        usually for mapping pair end data
        '''
        return {"pairs": [ [os.path.join(self.target_dir, self.id + "_treat_rep" + str(num+1)) + "pair1",
                  os.path.join(self.target_dir, self.id + "_treat_rep" + str(num+1)) + "pair2"]
                 for num in range(len(self.treatment_raws)) ],
                "reps": [os.path.join(self.target_dir,
                    self.id + "_treat_rep" + str(num+1)
                ) for num in range(len(self.treatment_raws))]}

    @property
    def control_pair_targets(self):
        '''pairs: for [[rep1_pair1, rep1_pair2]],
        usually for evaluating read quality
        reps: for [rep1, rep2],
        usually for mapping pair end data
        '''
        return {"pairs": [ [os.path.join(self.target_dir, self.id + "_control_rep" + str(num+1)) + "pair1",
                  os.path.join(self.target_dir, self.id + "_control_rep" + str(num+1)) + "pair2"]
                 for num in range(len(self.control_raws)) ],
                "reps": [os.path.join(self.target_dir,
                    self.id + "_control_rep" + str(num+1)
                ) for num in range(len(self.control_raws))]}

    @property
    def control_targets(self):
        if not self.pe:
            return self.control_single_targets
        else:
            return self.control_pair_targets["reps"]

    @property
    def control_single_targets(self):
        return [os.path.join(self.target_dir,
            self.id + "_control_rep" + str(num+1)
        ) for num in range(len(self.control_raws))]

    @property
    def sample_targets(self):
        return self.treatment_targets + self.control_targets

    def _base(self, path):
        return os.path.basename(path)

    @property
    def treatment_bases(self):
        return [self._base(i) for i in self.treatment_targets]

    @property
    def control_bases(self):
        return [self._base(i) for i in self.control_targets]

    @property
    def sample_bases(self):
        return [self._base(i) for i in self.sample_targets]

    def category(self, category_name):
        target_path = os.path.join(self.target_dir, category_name)
        if not os.path.exists(target_path):
            os.makedirs(target_path)
        return target_path

    @property
    def log(self):
        log_path = os.path.join(self.target_dir, 'log')
        if not os.path.exists(log_path):
            os.makedirs(log_path)
        logger = logging.getLogger();
        handler = logging.FileHandler(os.path.join(log_path, self.id + '.log'))
        logger.addHandler(handler);
        logger.setLevel(logging.NOTSET);
        return logger
