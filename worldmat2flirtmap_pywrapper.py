from nipype.interfaces.matlab import MatlabCommand
from nipype.interfaces.base import TraitedSpec, BaseInterface, BaseInterfaceInputSpec, File, traits
import os
from string import Template

class worldmat2flirtmapInputSpec(BaseInterfaceInputSpec):
    worldmat = File(exists=True, mandatory=True,
                        desc = 'matlab format matfile')
    src = File(exists=True, mandatory=True,
                    desc = 'source image')
    trg = File(exists=True, mandatory=True,
                    desc = 'target image')
    output_file = File(exists=False, mandatory=True,
                        desc = 'name of output file')


class worldmat2flirtmapOutputSpec(TraitedSpec):
    output_file = File(exists=True)

class worldmat2flirtmap(BaseInterface):
    input_spec = worldmat2flirtmapInputSpec
    output_spec = worldmat2flirtmapOutputSpec

    def _run_interface(self, runtime):
        d = dict(worldmat=self.inputs.worldmat,
        src=self.inputs.src,
        trg=self.inputs.trg,
        output_file=self.inputs.output_file)
        #this is your MATLAB code template
        script = Template("""worldmat = '$worldmat';
        src = '$src';
        trg = '$trg';
        output_file = '$output_file'
        worldmat2flirtmap(worldmat, src, trg, output_file);
        exit;
        """).substitute(d)

        # mfile = True  will create an .m file with your script and executed.
        # Alternatively
        # mfile can be set to False which will cause the matlab code to be
        # passed
        # as a commandline argument to the matlab executable
        # (without creating any files).
        # This, however, is less reliable and harder to debug
        # (code will be reduced to
        # a single line and stripped of any comments).

        mlab = MatlabCommand(script=script, mfile=True)
        result = mlab.run()
        return result.runtime

    def _list_outputs(self):
        outputs = self._outputs().get()
        outputs['output_file'] = os.path.abspath(self.inputs.output_file)
        return outputs
