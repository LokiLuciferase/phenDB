class PicaResultForUI:

    def __init__(self, job):
        self.job = job

    '''def __calc

def picaResultForUI(job):
    resultList = []
    all_bins = BinInJob.objects.filter(job=job)
    for bin_obj in all_bins:
        resultsList += PicaResult.objects.filter(bin=bin_obj.bin).select_related('bin') '''
