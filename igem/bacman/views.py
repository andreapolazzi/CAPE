from django.shortcuts import render

# Create your views here.
def homepage(request):
    return render(request, 'bacman/bacman-webpage.html')

def contact(request):
    return render(request, 'bacman/contactBacman.html')

def help(request):
    return render(request, 'bacman/helpBacman.html')