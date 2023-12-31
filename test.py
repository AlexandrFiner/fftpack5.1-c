import numpy as np

def test():
    a = '[1.000000 0.999980 0.999921 0.999822 0.999684 0.999507 0.999289 0.999033 0.998737 0.998402 0.998027 0.997612 0.997159 0.996666 0.996134 0.995562 0.994951 0.994301 0.993611 0.992883 0.992115 0.991308 0.990461 0.989576 0.988652 0.987688 0.986686 0.985645 0.984564 0.983445 0.982287 0.981090 0.979855 0.978581 0.977268 0.975917 0.974527 0.973099 0.971632 0.970127 0.968583 0.967001 0.965382 0.963724 0.962028 0.960294 0.958522 0.956712 0.954865 0.952979 0.951057 0.949096 0.947098 0.945063 0.942991 0.940881 0.938734 0.936550 0.934329 0.932071 0.929776 0.927445 0.925077 0.922673 0.920232 0.917755 0.915241 0.912692 0.910106 0.907484 0.904827 0.902134 0.899405 0.896641 0.893841 0.891007 0.888136 0.885231 0.882291 0.879316 0.876307 0.873262 0.870184 0.867071 0.863923 0.860742 0.857527 0.854277 0.850994 0.847678 0.844328 0.840945 0.837528 0.834078 0.830596 0.827081 0.823533 0.819952 0.816339 0.812694 0.809017 0.805308 0.801567 0.797794 0.793990 0.790155 0.786288 0.782391 0.778462 0.774503 0.770513 0.766493 0.762443 0.758362 0.754251 0.750111 0.745941 0.741742 0.737513 0.733255 0.728969 0.724653 0.720309 0.715936 0.711536 0.707107 0.702650 0.698165 0.693653 0.689114 0.684547 0.679953 0.675333 0.670686 0.666012 0.661312 0.656586 0.651834 0.647056 0.642253 0.637424 0.632570 0.627691 0.622788 0.617860 0.612907 0.607930 0.602930 0.597905 0.592857 0.587785 0.582690 0.577573 0.572432 0.567269 0.562083 0.556876 0.551646 0.546394 0.541121 0.535827 0.530511 0.525175 0.519817 0.514440 0.509041 0.503623 0.498185 0.492727 0.487250 0.481754 0.476238 0.470704 0.465151 0.459580 0.453991 0.448383 0.442758 0.437116 0.431456 0.425779 0.420086 0.414376 0.408649 0.402906 0.397148 0.391374 0.385584 0.379779 0.373959 0.368125 0.362275 0.356412 0.350534 0.344643 0.338738 0.332820 0.326888 0.320944 0.314986 0.309017 0.303035 0.297042 0.291036 0.285019 0.278991 0.272952 0.266902 0.260842 0.254771 0.248690 0.242599 0.236499 0.230389 0.224271 0.218143 0.212007 0.205863 0.199710 0.193549 0.187381 0.181206 0.175023 0.168833 0.162637 0.156434 0.150226 0.144011 0.137790 0.131564 0.125333 0.119097 0.112856 0.106611 0.100362 0.094108 0.087851 0.081591 0.075327 0.069060 0.062790 0.056519 0.050244 0.043968 0.037690 0.031411 0.025130 0.018848 0.012566 0.006283 1.000000 0.999921 0.999684 0.999289 0.998737 0.998027 0.997159 0.996134 0.994951 0.993611 0.992115 0.990461 0.988652 0.986686 0.984564 0.982287 0.979855 0.977268 0.974527 0.971632 0.968583 0.965382 0.962028 0.958522 0.954865 0.951057 0.947098 0.942991 0.938734 0.934329 0.929776 0.925077 0.920232 0.915241 0.910106 0.904827 0.899405 0.893841 0.888136 0.882291 0.876307 0.870184 0.863923 0.857527 0.850994 0.844328 0.837528 0.830596 0.823533 0.816339 0.809017 0.801567 0.793990 0.786288 0.778462 0.770513 0.762443 0.754251 0.745941 0.737513 0.728969 0.720309 0.711536 0.702650 0.693653 0.684547 0.675333 0.666012 0.656586 0.647056 0.637424 0.627691 0.617860 0.607930 0.597905 0.587785 0.577573 0.567269 0.556876 0.546394 0.535827 0.525175 0.514440 0.503623 0.492727 0.481754 0.470704 0.459580 0.448383 0.437116 0.425779 0.414376 0.402906 0.391374 0.379779 0.368125 0.356412 0.344643 0.332820 0.320944 0.309017 0.297042 0.285019 0.272952 0.260842 0.248690 0.236499 0.224271 0.212007 0.199710 0.187381 0.175023 0.162637 0.150226 0.137790 0.125333 0.112856 0.100362 0.087851 0.075327 0.062790 0.050244 0.037690 0.025130 0.012566 -0.000000 -0.012566 -0.025130 -0.037690 -0.050244 -0.062791 -0.075327 -0.087851 -0.100362 -0.112856 -0.125333 -0.137790 -0.150226 -0.162637 -0.175023 -0.187381 -0.199710 -0.212007 -0.224271 -0.236499 -0.248690 -0.260841 -0.272952 -0.285019 -0.297042 -0.309017 -0.320944 -0.332820 -0.344643 -0.356412 -0.368125 -0.379779 -0.391374 -0.402907 -0.414376 -0.425779 -0.437116 -0.448383 -0.459580 -0.470704 -0.481754 -0.492727 -0.503623 -0.514440 -0.525175 -0.535827 -0.546394 -0.556876 -0.567269 -0.577573 -0.587785 -0.597905 -0.607930 -0.617860 -0.627691 -0.637424 -0.647056 -0.656586 -0.666012 -0.675333 -0.684547 -0.693653 -0.702650 -0.711536 -0.720309 -0.728969 -0.737513 -0.745941 -0.754251 -0.762443 -0.770513 -0.778462 -0.786288 -0.793990 -0.801567 -0.809017 -0.816339 -0.823533 -0.830596 -0.837528 -0.844328 -0.850995 -0.857527 -0.863923 -0.870184 -0.876307 -0.882291 -0.888136 -0.893841 -0.899405 -0.904827 -0.910106 -0.915241 -0.920232 -0.925077 -0.929777 -0.934329 -0.938734 -0.942991 -0.947098 -0.951056 -0.954865 -0.958522 -0.962028 -0.965382 -0.968583 -0.971632 -0.974527 -0.977268 -0.979855 -0.982287 -0.984564 -0.986686 -0.988652 -0.990461 -0.992115 -0.993611 -0.994951 -0.996134 -0.997159 -0.998027 -0.998737 -0.999289 -0.999684 -0.999921 1.000000 0.999822 0.999289 0.998402 0.997159 0.995562 0.993611 0.991308 0.988652 0.985645 0.982287 0.978581 0.974527 0.970127 0.965382 0.960294 0.954865 0.949096 0.942991 0.936550 0.929776 0.922673 0.915241 0.907484 0.899405 0.891007 0.882291 0.873262 0.863923 0.854277 0.844328 0.834078 0.823533 0.812694 0.801567 0.790155 0.778462 0.766493 0.754251 0.741742 0.728969 0.715936 0.702650 0.689114 0.675333 0.661312 0.647056 0.632570 0.617860 0.602930 0.587785 0.572432 0.556876 0.541121 0.525175 0.509041 0.492727 0.476238 0.459580 0.442758 0.425779 0.408649 0.391374 0.373959 0.356412 0.338738 0.320944 0.303035 0.285019 0.266902 0.248690 0.230389 0.212007 0.193549 0.175023 0.156434 0.137790 0.119097 0.100362 0.081591 0.062791 0.043968 0.025130 0.006283 -0.012566 -0.031411 -0.050244 -0.069060 -0.087851 -0.106611 -0.125333 -0.144011 -0.162637 -0.181206 -0.199710 -0.218143 -0.236499 -0.254771 -0.272952 -0.291036 -0.309017 -0.326888 -0.344643 -0.362275 -0.379779 -0.397148 -0.414375 -0.431456 -0.448383 -0.465151 -0.481754 -0.498185 -0.514440 -0.530511 -0.546394 -0.562083 -0.577573 -0.592857 -0.607930 -0.622788 -0.637424 -0.651834 -0.666012 -0.679953 -0.693653 -0.707107 -0.720309 -0.733255 -0.745941 -0.758362 -0.770513 -0.782391 -0.793990 -0.805308 -0.816339 -0.827080 -0.837528 -0.847678 -0.857527 -0.867071 -0.876307 -0.885231 -0.893841 -0.902134 -0.910106 -0.917755 -0.925077 -0.932071 -0.938734 -0.945063 -0.951056 -0.956712 -0.962028 -0.967001 -0.971632 -0.975917 -0.979855 -0.983445 -0.986686 -0.989576 -0.992115 -0.994301 -0.996134 -0.997612 -0.998737 -0.999507 -0.999921 -0.999980 -0.999684 -0.999033 -0.998027 -0.996666 -0.994951 -0.992883 -0.990461 -0.987688 -0.984564 -0.981090 -0.977268 -0.973099 -0.968583 -0.963724 -0.958522 -0.952979 -0.947098 -0.940881 -0.934329 -0.927445 -0.920232 -0.912692 -0.904827 -0.896641 -0.888137 -0.879316 -0.870184 -0.860742 -0.850995 -0.840945 -0.830596 -0.819952 -0.809017 -0.797794 -0.786288 -0.774503 -0.762443 -0.750111 -0.737513 -0.724653 -0.711536 -0.698165 -0.684547 -0.670686 -0.656586 -0.642253 -0.627691 -0.612907 -0.597905 -0.582691 -0.567269 -0.551646 -0.535827 -0.519817 -0.503623 -0.487250 -0.470704 -0.453991 -0.437116 -0.420086 -0.402907 -0.385584 -0.368125 -0.350535 -0.332820 -0.314987 -0.297042 -0.278991 -0.260842 -0.242599 -0.224271 -0.205863 -0.187381 -0.168834 -0.150226 -0.131564 -0.112857 -0.094108 -0.075327 -0.056518 -0.037690 -0.018849 0.000000 0.006283 0.012566 0.018848 0.025130 0.031411 0.037690 0.043968 0.050244 0.056519 0.062791 0.069060 0.075327 0.081591 0.087851 0.094108 0.100362 0.106611 0.112856 0.119097 0.125333 0.131564 0.137790 0.144011 0.150226 0.156434 0.162637 0.168833 0.175023 0.181206 0.187381 0.193549 0.199710 0.205863 0.212007 0.218143 0.224271 0.230389 0.236499 0.242599 0.248690 0.254771 0.260842 0.266902 0.272952 0.278991 0.285019 0.291036 0.297042 0.303035 0.309017 0.314987 0.320944 0.326888 0.332820 0.338738 0.344643 0.350534 0.356412 0.362275 0.368125 0.373959 0.379779 0.385584 0.391374 0.397148 0.402906 0.408649 0.414376 0.420086 0.425779 0.431456 0.437116 0.442758 0.448383 0.453991 0.459580 0.465151 0.470704 0.476238 0.481754 0.487250 0.492727 0.498185 0.503623 0.509041 0.514440 0.519817 0.525175 0.530511 0.535827 0.541121 0.546394 0.551646 0.556876 0.562083 0.567269 0.572432 0.577573 0.582691 0.587785 0.592857 0.597905 0.602930 0.607930 0.612907 0.617860 0.622788 0.627691 0.632570 0.637424 0.642253 0.647056 0.651834 0.656586 0.661312 0.666012 0.670686 0.675333 0.679953 0.684547 0.689114 0.693653 0.698165 0.702650 0.707107 0.711536 0.715936 0.720309 0.724653 0.728969 0.733255 0.737513 0.741742 0.745941 0.750111 0.754251 0.758362 0.762443 0.766493 0.770513 0.774503 0.778462 0.782391 0.786288 0.790155 0.793990 0.797794 0.801567 0.805308 0.809017 0.812694 0.816339 0.819952 0.823533 0.827081 0.830596 0.834078 0.837528 0.840945 0.844328 0.847678 0.850994 0.854277 0.857527 0.860742 0.863923 0.867071 0.870184 0.873262 0.876307 0.879316 0.882291 0.885231 0.888136 0.891007 0.893841 0.896641 0.899405 0.902134 0.904827 0.907484 0.910106 0.912692 0.915241 0.917755 0.920232 0.922673 0.925077 0.927445 0.929776 0.932071 0.934329 0.936550 0.938734 0.940881 0.942991 0.945063 0.947098 0.949096 0.951057 0.952979 0.954865 0.956712 0.958522 0.960294 0.962028 0.963724 0.965382 0.967001 0.968583 0.970127 0.971632 0.973099 0.974527 0.975917 0.977268 0.978581 0.979855 0.981090 0.982287 0.983445 0.984564 0.985645 0.986686 0.987688 0.988652 0.989576 0.990461 0.991308 0.992115 0.992883 0.993611 0.994301 0.994951 0.995562 0.996134 0.996666 0.997159 0.997613 0.998027 0.998402 0.998737 0.999033 0.999289 0.999507 0.999684 0.999822 0.999921 0.999980 0.000000 0.012566 0.025130 0.037690 0.050244 0.062791 0.075327 0.087851 0.100362 0.112856 0.125333 0.137790 0.150226 0.162637 0.175023 0.187381 0.199710 0.212007 0.224271 0.236499 0.248690 0.260842 0.272952 0.285019 0.297042 0.309017 0.320944 0.332820 0.344643 0.356412 0.368125 0.379779 0.391374 0.402906 0.414376 0.425779 0.437116 0.448383 0.459580 0.470704 0.481754 0.492727 0.503623 0.514440 0.525175 0.535827 0.546394 0.556876 0.567269 0.577573 0.587785 0.597905 0.607930 0.617860 0.627691 0.637424 0.647056 0.656586 0.666012 0.675333 0.684547 0.693653 0.702650 0.711536 0.720309 0.728969 0.737513 0.745941 0.754251 0.762443 0.770513 0.778462 0.786288 0.793990 0.801567 0.809017 0.816339 0.823533 0.830596 0.837528 0.844328 0.850994 0.857527 0.863923 0.870184 0.876307 0.882291 0.888136 0.893841 0.899405 0.904827 0.910106 0.915241 0.920232 0.925077 0.929776 0.934329 0.938734 0.942991 0.947098 0.951057 0.954865 0.958522 0.962028 0.965382 0.968583 0.971632 0.974527 0.977268 0.979855 0.982287 0.984564 0.986686 0.988652 0.990461 0.992115 0.993611 0.994951 0.996134 0.997159 0.998027 0.998737 0.999289 0.999684 0.999921 1.000000 0.999921 0.999684 0.999289 0.998737 0.998027 0.997159 0.996134 0.994951 0.993611 0.992115 0.990461 0.988652 0.986686 0.984564 0.982287 0.979855 0.977268 0.974527 0.971632 0.968583 0.965382 0.962028 0.958522 0.954865 0.951056 0.947098 0.942990 0.938734 0.934329 0.929776 0.925077 0.920232 0.915241 0.910106 0.904827 0.899405 0.893841 0.888136 0.882291 0.876307 0.870184 0.863923 0.857527 0.850995 0.844328 0.837528 0.830596 0.823533 0.816339 0.809017 0.801567 0.793990 0.786288 0.778462 0.770513 0.762442 0.754251 0.745941 0.737513 0.728969 0.720309 0.711536 0.702650 0.693653 0.684547 0.675333 0.666012 0.656586 0.647056 0.637424 0.627691 0.617860 0.607930 0.597905 0.587785 0.577573 0.567269 0.556876 0.546394 0.535827 0.525175 0.514440 0.503623 0.492727 0.481754 0.470704 0.459580 0.448383 0.437116 0.425779 0.414375 0.402906 0.391374 0.379779 0.368124 0.356412 0.344643 0.332819 0.320944 0.309017 0.297041 0.285019 0.272952 0.260841 0.248690 0.236499 0.224271 0.212007 0.199710 0.187381 0.175023 0.162637 0.150225 0.137790 0.125333 0.112856 0.100362 0.087851 0.075327 0.062790 0.050244 0.037690 0.025130 0.012566 0.000000 0.018848 0.037690 0.056519 0.075327 0.094108 0.112856 0.131564 0.150226 0.168833 0.187381 0.205863 0.224271 0.242599 0.260842 0.278991 0.297042 0.314986 0.332820 0.350534 0.368125 0.385584 0.402906 0.420086 0.437116 0.453990 0.470704 0.487250 0.503623 0.519817 0.535827 0.551646 0.567269 0.582690 0.597905 0.612907 0.627691 0.642253 0.656586 0.670686 0.684547 0.698165 0.711536 0.724653 0.737513 0.750111 0.762443 0.774503 0.786288 0.797794 0.809017 0.819952 0.830596 0.840945 0.850994 0.860742 0.870184 0.879316 0.888136 0.896641 0.904827 0.912692 0.920232 0.927445 0.934329 0.940881 0.947098 0.952979 0.958522 0.963724 0.968583 0.973099 0.977268 0.981090 0.984564 0.987688 0.990461 0.992883 0.994951 0.996666 0.998027 0.999033 0.999684 0.999980 0.999921 0.999507 0.998737 0.997613 0.996134 0.994301 0.992115 0.989576 0.986686 0.983445 0.979855 0.975917 0.971632 0.967001 0.962028 0.956712 0.951057 0.945063 0.938734 0.932071 0.925077 0.917755 0.910106 0.902134 0.893841 0.885231 0.876307 0.867071 0.857527 0.847678 0.837528 0.827081 0.816339 0.805308 0.793990 0.782391 0.770513 0.758362 0.745941 0.733255 0.720309 0.707107 0.693653 0.679953 0.666012 0.651834 0.637424 0.622788 0.607930 0.592857 0.577573 0.562083 0.546394 0.530511 0.514440 0.498185 0.481754 0.465151 0.448383 0.431456 0.414376 0.397148 0.379779 0.362275 0.344643 0.326888 0.309017 0.291036 0.272952 0.254771 0.236499 0.218143 0.199710 0.181206 0.162637 0.144011 0.125333 0.106611 0.087851 0.069060 0.050244 0.031411 0.012566 -0.006283 -0.025130 -0.043968 -0.062790 -0.081591 -0.100362 -0.119097 -0.137790 -0.156434 -0.175023 -0.193549 -0.212007 -0.230389 -0.248690 -0.266902 -0.285019 -0.303035 -0.320944 -0.338738 -0.356412 -0.373959 -0.391374 -0.408649 -0.425779 -0.442758 -0.459580 -0.476238 -0.492727 -0.509041 -0.525175 -0.541121 -0.556876 -0.572432 -0.587785 -0.602929 -0.617860 -0.632570 -0.647056 -0.661312 -0.675333 -0.689114 -0.702650 -0.715936 -0.728969 -0.741742 -0.754251 -0.766493 -0.778462 -0.790155 -0.801567 -0.812694 -0.823532 -0.834078 -0.844328 -0.854277 -0.863923 -0.873262 -0.882291 -0.891006 -0.899405 -0.907484 -0.915241 -0.922673 -0.929776 -0.936550 -0.942991 -0.949096 -0.954865 -0.960294 -0.965382 -0.970127 -0.974527 -0.978581 -0.982287 -0.985645 -0.988652 -0.991308 -0.993611 -0.995562 -0.997159 -0.998402 -0.999289 -0.999822 1.000000 0.999684 0.998737 0.997159 0.994951 0.992115 0.988652 0.984564 0.979855 0.974527 0.968583 0.962028 0.954865 0.947098 0.938734 0.929776 0.920232 0.910106 0.899405 0.888136 0.876307 0.863923 0.850994 0.837528 0.823533 0.809017 0.793990 0.778462 0.762443 0.745941 0.728969 0.711536 0.693653 0.675333 0.656586 0.637424 0.617860 0.597905 0.577573 0.556876 0.535827 0.514440 0.492727 0.470704 0.448383 0.425779 0.402906 0.379779 0.356412 0.332820 0.309017 0.285019 0.260842 0.236499 0.212007 0.187381 0.162637 0.137790 0.112856 0.087851 0.062790 0.037690 0.012566 -0.012566 -0.037690 -0.062791 -0.087851 -0.112856 -0.137790 -0.162637 -0.187381 -0.212007 -0.236499 -0.260841 -0.285019 -0.309017 -0.332820 -0.356412 -0.379779 -0.402907 -0.425779 -0.448383 -0.470704 -0.492727 -0.514440 -0.535827 -0.556876 -0.577573 -0.597905 -0.617860 -0.637424 -0.656586 -0.675333 -0.693653 -0.711536 -0.728969 -0.745941 -0.762443 -0.778462 -0.793990 -0.809017 -0.823533 -0.837528 -0.850995 -0.863923 -0.876307 -0.888136 -0.899405 -0.910106 -0.920232 -0.929777 -0.938734 -0.947098 -0.954865 -0.962028 -0.968583 -0.974527 -0.979855 -0.984564 -0.988652 -0.992115 -0.994951 -0.997159 -0.998737 -0.999684 0.000000 0.025130 0.050244 0.075327 0.100362 0.125333 0.150226 0.175023 0.199710 0.224271 0.248690 0.272952 0.297042 0.320944 0.344643 0.368125 0.391374 0.414376 0.437116 0.459580 0.481754 0.503623 0.525175 0.546394 0.567269 0.587785 0.607930 0.627691 0.647056 0.666012 0.684547 0.702650 0.720309 0.737513 0.754251 0.770513 0.786288 0.801567 0.816339 0.830596 0.844328 0.857527 0.870184 0.882291 0.893841 0.904827 0.915241 0.925077 0.934329 0.942991 0.951057 0.958522 0.965382 0.971632 0.977268 0.982287 0.986686 0.990461 0.993611 0.996134 0.998027 0.999289 0.999921 0.999921 0.999289 0.998027 0.996134 0.993611 0.990461 0.986686 0.982287 0.977268 0.971632 0.965382 0.958522 0.951056 0.942990 0.934329 0.925077 0.915241 0.904827 0.893841 0.882291 0.870184 0.857527 0.844328 0.830596 0.816339 0.801567 0.786288 0.770513 0.754251 0.737513 0.720309 0.702650 0.684547 0.666012 0.647056 0.627691 0.607930 0.587785 0.567269 0.546394 0.525175 0.503623 0.481754 0.459580 0.437116 0.414375 0.391374 0.368124 0.344643 0.320944 0.297041 0.272952 0.248690 0.224271 0.199710 0.175023 0.150225 0.125333 0.100362 0.075327 0.050244 0.025130 1.000000 0.998737 0.994951 0.988652 0.979855 0.968583 0.954865 0.938734 0.920232 0.899405 0.876307 0.850994 0.823533 0.793990 0.762443 0.728969 0.693653 0.656586 0.617860 0.577573 0.535827 0.492727 0.448383 0.402906 0.356412 1.000000 0.994951 0.979855 0.954865 0.920232 0.876307 0.823533 0.762443 0.693653 0.617860 0.535827 0.448383 0.356412 0.260842 0.162637 0.062790 -0.037690 -0.137790 -0.236499 -0.332820 -0.425779 -0.514440 -0.597905 -0.675333 -0.745941 1.000000 0.988652 0.954865 0.899405 0.823533 0.728969 0.617860 0.492727 0.356412 0.212007 0.062791 -0.087851 -0.236499 -0.379779 -0.514440 -0.637424 -0.745941 -0.837528 -0.910106 -0.962028 -0.992115 -0.999684 -0.984564 -0.947098 -0.888137 1.000000 0.979855 0.920232 0.823533 0.693653 0.535827 0.356412 0.162637 -0.037690 -0.236499 -0.425779 -0.597905 -0.745941 -0.863923 -0.947098 -0.992115 -0.997159 -0.962028 -0.888136 -0.778462 -0.637424 -0.470704 -0.285019 -0.087851 0.112856 0.000000 0.050244 0.100362 0.150226 0.199710 0.248690 0.297042 0.344643 0.391374 0.437116 0.481754 0.525175 0.567269 0.607930 0.647056 0.684547 0.720309 0.754251 0.786288 0.816339 0.844328 0.870184 0.893841 0.915241 0.934329 0.000000 0.100362 0.199710 0.297042 0.391374 0.481754 0.567269 0.647056 0.720309 0.786288 0.844328 0.893841 0.934329 0.965382 0.986686 0.998027 0.999289 0.990461 0.971632 0.942990 0.904827 0.857527 0.801567 0.737513 0.666012 0.000000 0.150226 0.297042 0.437116 0.567269 0.684547 0.786288 0.870184 0.934329 0.977268 0.998027 0.996134 0.971632 0.925077 0.857527 0.770513 0.666012 0.546394 0.414376 0.272952 0.125333 -0.025130 -0.175023 -0.320944 -0.459580 0.000000 0.199710 0.391374 0.567269 0.720309 0.844328 0.934329 0.986686 0.999289 0.971632 0.904827 0.801567 0.666012 0.503623 0.320944 0.125333 -0.075327 -0.272952 -0.459580 -0.627692 -0.770513 -0.882291 -0.958522 -0.996134 -0.993611 1.000000 0.968583 0.876307 0.728969 0.535827 1.000000 0.876307 0.535827 0.062790 -0.425779 1.000000 0.728969 0.062790 -0.637424 -0.992115 1.000000 0.535827 -0.425779 -0.992115 -0.637424 0.000000 0.248690 0.481754 0.684547 0.844328 0.000000 0.481754 0.844328 0.998027 0.904827 0.000000 0.684547 0.998027 0.770513 0.125333 0.000000 0.844328 0.904827 0.125333 -0.770513 1.000000 1.000000 1.000000 1.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 5.000000 4.000000 2.000000 5.000000 5.000000 5.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 ]'

def CFFT1I(signal):
    # Обратное преобразование Фурье
    result = np.fft.ifft(signal)
    return result

def CFFT1B(signal):
    # Прямое преобразование Фурье
    spectrum = np.fft.fft(signal)
    return spectrum

def CFFT1F(signal):
    # Прямое преобразование Фурье
    result = np.fft.fft(signal)
    return result

# Пример использования
if __name__ == "__main__":
    # Генерация входных данных
    signal = np.array([1, 2, 3, 4])

    # Прямое преобразование Фурье
    result_fft = CFFT1F(signal)
    print("Input signal:", signal)
    print("FFT result:", result_fft)

    # Обратное преобразование Фурье
    result_ifft = CFFT1I(result_fft)
    print("IFFT result:", result_ifft)

    # Прямое преобразование Фурье и обратное преобразование Фурье
    result_fft_ifft = CFFT1B(signal)
    print("FFT result:", result_fft_ifft)
    print("IFFT result:", CFFT1I(result_fft_ifft))