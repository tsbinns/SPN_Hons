from neuron import h
h.nrnmpi_init()

pc = h.ParallelContext()
print ("I am %d of %d" % (pc.id(), pc.nhost()))

pc.barrier()
h.quit()