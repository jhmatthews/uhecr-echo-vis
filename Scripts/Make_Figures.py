import fig1, fig2, fig3, fig4, fig5, argparse

def MakeFigures(load=True):
    '''
    make all the figures
    '''
    fig1.run()
    fig2.run()
    fig3.run(load=load)
    fig4.run(load=load)
    fig5.run()

def parse_args():
    '''parse arguments'''
    parser = argparse.ArgumentParser(
        prog='MakeFigures',
        description='make figures 1-5 for Taylor et al. UHECR echoes paper')
    parser.add_argument('-r', '--raw', action='store_true',
                        default=False, help="Make movie frames")
    args = parser.parse_args()
    return args

if __name__ == "__main__":
    args = parse_args()
    load = (args.raw == False)
    MakeFigures(load=load)