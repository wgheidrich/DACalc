import io
import sys
import re
from contextlib import redirect_stderr
from contextlib import redirect_stdout
from ipykernel.kernelbase import Kernel
from dacalc.calculator import parse

import numpy
import matplotlib.pyplot as plt
import base64
import io


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
        
        # although the parser can handle multi-line input, we go line
        # by line so that error messages are in the right spot
        results = []
        for line in code.splitlines():
            with io.StringIO() as buf, redirect_stdout(buf), \
                 io.StringIO() as err_buf, redirect_stderr(err_buf):
                try:
                    results += parse(line)
                except Exception as exc:
                    print(exc, file=sys.stderr)
                    error = True
                error_msg = err_buf.getvalue()
                if error_msg != "":
                    output = "<font color=#FF5733>" + error_msg + "</font>"
                data = buf.getvalue()
                if isinstance(data,str):
                    output += data

        if isinstance(output,str) and not silent:
            # output is HTML so that special characters and tags can
            # be used in user messages
            output = "<tt>" + re.sub("\n", "<br>", output) + "</tt>"
            self.send_response(self.iopub_socket,
                               'execute_result',
                               {"execution_count": self.execution_count,
                                "data": {"text/html": output},
                                "metadata": {}})

        # look for images and paged text in the results and show them
        help_msg = ""
        for item in results:
            if isinstance(item, str):
                # concatenate all help messages into a single string
                if item != "":
                    help_msg += "\n\n---\n\n" + item
            elif isinstance(item, numpy.ndarray):
                # convert arrays to images and send display message
                buf = io.BytesIO()
                plt.imsave(buf, arr=item, format='png')
                buf.seek(0)
                png_string = base64.b64encode(buf.read()).decode('utf-8')
                self.send_response(self.iopub_socket, 'display_data',
                                   {"data": {"image/png": png_string},
                                    "metadata": {}})
            
        if help_msg != "":
            popup = [{"source": "page", "data": {"text/plain": help_msg}}]
        else:
            popup = []

        return {'status': 'error' if error else 'ok',
                # The base class increments the execution count
                'execution_count': self.execution_count,
                'payload': popup,
                'user_expressions': {}}
