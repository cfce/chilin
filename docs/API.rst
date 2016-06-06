chilin2
=======

extend your pipeline
---------------------
ChiLin is written with python and shell, it's easy to extend to fit your need.

* in chilin2/modules/, `mkdir chilin2/modules/my_pipe && touch chilin2/modules/my_pipe/__init__.py`
* write your command into chilin2/modules/my_pipe/__init__.py like following:

.. code-block:: python

    from samflow.command import ShellCommand
    from samflow.workflow import attach_back
    def tail(workflow, conf):
        tail = attach_back(workflow,
        ShellCommand(
        """{tool} -50 {input[peaks]} > {output[p]}
        """,
        tool = "sort",
        input = {"peaks": conf.prefix + "_sort_peaks.narrowPeak"},
        output = {"p":conf.prefix + "_sort_peaks_tail_50.narrowPeak"},
        name = "get tail lines"))


.. code-block:: python

    from chilin2.modules.my_pipe import tail
    ....
    #at line 349
    if need_run(13):
       bld.build(tail)


.. :mod:`chilin2` -- utils package
.. ---------------------------------
.. .. automodule:: chilin2
.. 
.. .. automodule:: chilin2.modules
.. 
.. .. automodule:: chilin2.modules.bwa
.. 
.. .. automodule:: chilin2.modules.config
.. 
.. .. automodule:: chilin2.modules.config.config
.. 
.. .. automodule:: chilin2.modules.config.helpers
.. 
.. .. automodule:: chilin2.modules.fastqc.dc
.. 
.. .. automodule:: chilin2.modules.fastqc.qc
.. 
.. .. automodule:: chilin2.modules.fastqc.tex
.. 
.. .. automodule:: chilin2.modules.bwa.dc
.. 
.. .. automodule:: chilin2.modules.bwa.qc
.. 
.. .. automethod:: chilin2.modules.bwa.dc.bwa
.. 
.. .. automodule:: chilin2.modules.macs
.. 
.. .. automodule:: chilin2.modules.macs.dc
.. 
.. .. automethod:: chilin2.modules.macs.dc.macs2
.. 
.. .. automodule:: chilin2.modules.ceas
.. 
.. .. automodule:: chilin2.modules.ceas.dc
.. 
.. .. automodule:: chilin2.modules.ceas.qc
.. 
.. .. automodule:: chilin2.modules.frip
.. 
.. .. automodule:: chilin2.modules.frip.dc
.. 
.. .. automodule:: chilin2.modules.frip.qc
.. 
.. .. automodule:: chilin2.modules.library
.. 
.. .. automodule:: chilin2.modules.library.dc
.. 
.. .. automodule:: chilin2.modules.library.qc
.. 
.. .. automodule:: chilin2.modules.enrichment
.. 
.. .. automodule:: chilin2.modules.enrichment.dc
.. 
.. .. automodule:: chilin2.modules.macs2_fragment
.. 
.. .. automodule:: chilin2.modules.macs2_fragment.dc
.. 
.. .. automodule:: chilin2.modules.macs2_fragment.qc
.. 
.. .. automodule:: chilin2.modules.replicates
.. 
.. .. automodule:: chilin2.modules.replicates.dc
.. 
.. .. automodule:: chilin2.modules.replicates.qc
.. 
.. .. automodule:: chilin2.modules.conservation
.. 
.. .. automodule:: chilin2.modules.conservation.dc
.. 
.. .. automodule:: chilin2.modules.conservation.qc
.. 
.. .. automodule:: chilin2.modules.contamination
.. 
.. .. automodule:: chilin2.modules.contamination.dc
.. 
.. .. automodule:: chilin2.modules.contamination.qc
.. 
.. .. automodule:: chilin2.modules.contamination.tex
.. 
.. .. automodule:: chilin2.modules.regulatory
.. 
.. .. automodule:: chilin2.modules.regulatory.dc
.. 
.. .. automodule:: chilin2.modules.mdseqpos
.. 
.. .. automodule:: chilin2.modules.mdseqpos.dc
.. 
.. .. automodule:: chilin2.modules.mdseqpos.qc
.. 
.. .. automodule:: chilin2.modules.mdseqpos.tex
.. 
.. .. automodule:: chilin2.modules.summary
.. 
.. .. automodule:: chilin2.modules.summary.qc_json_summary
.. 
.. .. automodule:: chilin2.modules.summary.qc_summary_table
.. 
.. :mod:`samflow` -- utils package
.. ---------------------------------
.. 
.. .. automodule:: samflow
.. 
.. .. autoclass:: samflow.command.ShellCommand
.. 
.. .. autoclass:: samflow.workflow.Workflow
.. 
.. .. autoclass:: samflow.command.PythonCommand
