import io
import sys
from contextlib import redirect_stderr
from contextlib import redirect_stdout
from ipykernel.kernelbase import Kernel
from dacalc.calculator import parse

class DAKernel(Kernel):
    implementation = 'DimensionalAnalysis'
    implementation_version = '1.0'
    language = 'no-op'
    language_version = '0.1'
    language_info = {
        'name': 'dacalc',
        'mimetype': 'text/plain',
        'file_extension': '.da',
    }
    banner = '''Dimensional Analysis Calculator

try '?' for help...
'''

        
    
    def do_execute(self, code, silent,
                   store_history=True, user_expressions=None,
                   allow_stdin=False):
        error = False
        output = ""
        help_msg = ""
        for line in code.splitlines():
            with io.StringIO() as buf, redirect_stdout(buf), io.StringIO() as err_buf, redirect_stderr(err_buf):
                try:
                    help_msg += parse(line)
                except Exception as exc:
                    print(exc,file = sys.stderr)
                    error = True
                error_msg = err_buf.getvalue()
                if error_msg != "":
                    output += '\u001b[38;5;160m'+error_msg+"\u001b[0m"
                output += buf.getvalue()
        
        if not silent:
            self.send_response(self.iopub_socket,
                               'execute_result',
                               {"execution_count": self.execution_count,
                                "data": {"text/plain": output},
                                "metadata": {}})

        if type(help_msg) != str or help_msg == "":
            print("no help message")
            popup = []
        else:
            print("trying help message popup")
            popup = [{"source" : "page",
                       "data" : {"text/plain" : help_msg}
            }]
                
        return {'status': 'error' if error else 'ok',
                # The base class increments the execution count
                'execution_count': self.execution_count,
                'payload': popup,
                'user_expressions': {},
               }
