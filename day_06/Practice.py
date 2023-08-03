# -*- coding: utf-8 -*-
"""
Created on Tue Aug  1 01:34:31 2023

@author: USER
"""
#%%
sequence = [1, 2, 3]
sequence*3  # ? what does this do
#%%
negative_ned = (False, False, False, False, True, False)  # He only likes to be negative
for elem in negative_ned:
    if not elem:
        print("Ned said no again")
    else:
        print("Ned did the impossible.")
        
#%%
# range() makes an integer iterable
simple_comprehension = [x for x in range(4)]  # [0, 1, 2, 3]
simple_comprehension
#%%
comprehended = [2*x for x in range(10) if x**2 > 3]  # limit to 10
comprehended