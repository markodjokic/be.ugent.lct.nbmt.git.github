import static org.junit.Assert.*;

import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;


public class SimpleTest {

	@BeforeClass
	public static void setUpBeforeClass() throws Exception {
	}

	@AfterClass
	public static void tearDownAfterClass() throws Exception {
	}

	@Before
	public void setUp() throws Exception {
	}

	@After
	public void tearDown() throws Exception {
	}

	@Test
	public void testSum() {
		int[] vec = {4,3,2};
		int sum = Simple.sum(vec);
		assertTrue(sum==9);
	}

	@Test
	public void testProd() {
		int[] vec = {4,3,2};
		int sum = Simple.prod(vec);
		assertTrue(sum==24);
	}

}
