#include <stdio.h>
#include <string.h>

/****************************************************************************************/

int main(void)
{
	FILE *fromFile;
	FILE *toFile;
	int x, y;
	char buffer[1000], fromFileName[30], toFileName[30], molType[10], fileNamePrefix[20], numBuff1[5], numBuff2[5];

	strcpy(molType, "l");
	while(strcmp(molType, "") != 0)
	{
		for(x = 1; x <= 190; x++)
		{
			for(y = 1; y < 10; y++)
			{
				strcpy(fileNamePrefix, molType);
				strcat(fileNamePrefix, "-active.mol2_");
				if(x < 10)
				strcat(fileNamePrefix, "0");
				strcpy(numBuff1, "");
				sprintf(numBuff1, "%d", x);
				strcat(fileNamePrefix, numBuff1);
				strcat(fileNamePrefix, "-");
				strcpy(numBuff2, "");
				sprintf(numBuff2, "%d", y);
				strcat(fileNamePrefix, numBuff2);
				strcpy(fromFileName, fileNamePrefix);
				strcat(fromFileName, ".sdf");
				strcpy(toFileName, fileNamePrefix);
				strcat(toFileName, ".mol2");

				fromFile = fopen(fromFileName, "r");
				if(fromFile != NULL)
					toFile = fopen(toFileName, "a");
				else
					toFile = NULL;

				if(fromFile != NULL && toFile != NULL)
				{
					printf("%s\t%s\n", fromFileName, toFileName);
					fscanf(fromFile, " %[^\n]", buffer);
					while(!feof(fromFile) && strstr(buffer, "END") == NULL)
						fscanf(fromFile, " %[^\n]", buffer);
					while(!feof(fromFile))
					{
						fprintf(toFile, "%s\n", buffer);
						fscanf(fromFile, " %[^\n]", buffer);
					}
				}
				if(fromFile != NULL)
					fclose(fromFile);
				if(toFile != NULL)
					fclose(toFile);
			}
		}
		if(strcmp(molType, "l") == 0)
			strcpy(molType, "r");
		else
			strcpy(molType, "");
	}
	

	printf("Done.\n");
	return 0;
} /* main */

/****************************************************************************************/
/****************************************************************************************/
/****************************************************************************************/
/****************************************************************************************/
/****************************************************************************************/
/****************************************************************************************/
/****************************************************************************************/
/****************************************************************************************/