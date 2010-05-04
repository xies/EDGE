/**
 * Fit a straight line to the data, of the form z = mX(x - bX), z = mY(y - bY).
 * These form the equation to a line, z = mX(x - bX), z = mY(y - bY),
 * where m is the slope and b is the x-intercept
 * This is a helper class that holds four public final fields: mX, bX, mY, and bY.
 * Note: this routine ignores NaN rows.
 */
public class StraightLineFit {
	public final float mX;
	public final float bX;
	public final float mY;
	public final float bY;
	
	public StraightLineFit(float[][] yx, float[] z) {
		if (yx.length != z.length) {
			System.err.println("Error in StraightLineFit constructor: xy and z arrays must be the same length.");
			mX = Float.NaN;
			bX = Float.NaN;
			mY = Float.NaN;
			bY = Float.NaN;
			return;
		}
		if (yx[0].length != 2) {
			System.err.println("Cannot create StraightLineFit - input must have dimensions n x 2");
			mX = Float.NaN;
			bX = Float.NaN;
			mY = Float.NaN;
			bY = Float.NaN;
			return;
		}
		int n = z.length;

		float[] x = new float[n];
		float[] y = new float[n];
		for (int i = 0; i < n; i++) {
			y[i] = yx[i][0];
			x[i] = yx[i][1];
		}
		
		// remove NaN rows
		int notNaN = 0;
		for (int i = 0; i < n; i++) 
			if (!Float.isNaN(x[i])) notNaN++;
			
		float[] x2 = new float[notNaN];
		float[] y2 = new float[notNaN];
		float[] z2 = new float[notNaN];
		
		notNaN = 0;
		for (int i = 0; i < n; i++) { 
			if (!Float.isNaN(x[i])) {
				x2[notNaN] = x[i];
				y2[notNaN] = y[i];
				z2[notNaN] = z[i];
				notNaN++;
			}
		}
		
		// finished removing NaN rows

		float[] results;
		results = lineFit(x2, z2);
		mX = results[0];
		bX = results[1];
		results = lineFit(y2, z2);
		mY = results[0];
		bY = results[1];
	}

	// formula from Wikipedia: Linear Regression (subsection 2.1.1: univariate linear case)
	// performs the fitting with switching x and z, to avoid near infinite slopes
	private float[] lineFit(float[] x, float[] z) {
		int n = x.length;
		float[] temp = x;
		x = z;
		z = temp;
		
		float sumX = 0;
		float sumZ = 0;
		float sumSqX = 0;
		float sumXZ = 0;
		for (int i = 0; i < n; i++) {
			sumX += x[i];
			sumZ += z[i];
			sumSqX += x[i] * x[i];
			sumXZ += x[i] * z[i];
		}

		float mResult; 
		float bResult;

		mResult = (n * sumXZ - sumX * sumZ) / (n * sumSqX - sumX * sumX);
		bResult = (sumZ - mResult * sumX) / n;
		
		if (mResult == 0.0) {
			mResult = Float.MAX_VALUE;
			bResult = z[0];
		}
		else {
			mResult = 1.0f / mResult;
		}
		
		float[] results = new float[2];
		results[0] = mResult;
		results[1] = bResult;
		return results;
	}

	public String toString() {
		return "StraightLineFit, mX = " + mX + ", bX = " + bX + ", mY = " + mY + ", bY = " + bY + ".";
	}
	
}

