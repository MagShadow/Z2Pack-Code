from datetime import datetime

settings = {'num_lines': 51,
            'pos_tol':  1e-2,
            'gap_tol': 0.05,
            'move_tol': 0.2,
            'iterator': range(50, 81, 4),
            'min_neighbour_dist': 1e-4,
            }
settings_strict = {'num_lines': 81,
                   'pos_tol':  1e-3,
                   'gap_tol': 0.001,
                   'move_tol': 0.4,
                   'iterator': range(80, 201, 5),
                   'min_neighbour_dist': 5e-6,
                   }
def nt():
    '''
    Format the time.
    '''
    return datetime.now().strftime("%y-%m-%d-%H-%M-%S")