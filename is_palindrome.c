#include <stdio.h>
#include <string.h>
#include <ctype.h>

//(clear && gcc -g is_palindrome.c -o is_palindrome -lm -D _GNU_SOURCE) && (./is_palindrome)

int isPalindrome(const char *word) {
    int i, j;

    for (i = 0, j = strlen(word) - 1; i < j; i++, j--) {
        // compare caracters regardless of letter case

        if (tolower(word[i]) != tolower(word[j])) {
            return 0;  // break as the word is not a palindrome
        }
    }

    return 1;  // word is palindrome
}

int main() {
    char word[100];

    printf("Enter a word : ");
    scanf("%s", word);

    if (isPalindrome(word)) {
        printf("%s is a palindrome.\n", word);
    } else {
        printf("%s is not palindrome.\n", word);
    }

    return 0;
}
