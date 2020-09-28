from ipykernel.kernelapp import IPKernelApp
from .kernel import DAKernel

IPKernelApp.launch_instance(kernel_class=DAKernel)
