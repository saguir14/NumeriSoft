from django import forms

class CreateNewList(forms.Form):
    x0 = forms.FloatField(label = "x0")
    delta = forms.FloatField(label = "delta")
    iterations = forms.IntegerField(label = "iter")
    func = forms.CharField(label = "function", max_length = 200 )
    #check = forms.BooleanField()

