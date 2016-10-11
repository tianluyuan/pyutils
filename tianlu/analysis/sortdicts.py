""" Module containing dicts for sorting data
"""
MCP_DICT_OLD = {'RUN1':'2010-02-water/magnet/beama',
                'RUN2water':'2010-11-water/magnet/beamb',
                'RUN2air':'2010-11-air/magnet/beamb',
                'RUN3':'2010-11-air/magnet/beamc',
                'RUN4water':'2010-11-water/magnet/beamc',
                'RUN4air':'2010-11-air/magnet/beamc',
                'RUN3sand':'2010-11-air/sand/beamc',
                'RUN4sand':'2010-11-water/sand/beamc'}

MCP_DICT_NEW = {'RUN1':'2010-02-water/magnet/run1',
                'RUN2water':'2010-11-water/magnet/run2',
                'RUN2air':'2010-11-air/magnet/run2',
                'RUN3':'2010-11-air/magnet/run3',
                'RUN4water':'2010-11-water/magnet/run4',
                'RUN4air':'2010-11-air/magnet/run4',
                'RUN3sand':'2010-11-air/sand/run3',
                'RUN4sand':'2010-11-water/sand/run4'}

RDP_DICT = {'RUN1':'4165-5115',
            'RUN2water':'6462-7663',
            'RUN2air':'7665-7754',
            'RUN3b':'8309-8453',
            'RUN3c':'8550-8753',
            'RUN4water':'8983-9413',
            'RUN4air':'9426-9796'}

P0DGEOM_DICT = {'RUN1':'2010-02-water',
                'RUN2water':'2010-11-water',
                'RUN2air':'2010-11-air',
                'RUN3':'2010-11-air',
                'RUN3b':'2010-11-air',
                'RUN3c':'2010-11-air',
                'RUN4water':'2010-11-water',
                'RUN4air':'2010-11-air'}

def rdp_to_mcp(run):
    """ converts a rdp run to a mcp run (only really affects run3)
    """
    conversion = {'RUN3b':'RUN3', 'RUN3c':'RUN3'}
    return conversion[run] if run in conversion else run

