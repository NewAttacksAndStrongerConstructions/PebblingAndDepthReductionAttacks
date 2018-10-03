using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ConsoleApplication1
{
    // Generic heap to help make things nice and fast
    class Heap<T> where T : IComparable
    {
        private T[] arr;
        private int size;
        private int numContents;

        public Heap(int start_size = 10)
        {
            if(start_size < 2)
                throw new ArgumentException("Starting array size must be at least 2");

            arr = new T[start_size];
            size = start_size;
            numContents = 0;
        }

        public T minValue()
        {
            if (numContents < 1)
                throw new IndexOutOfRangeException("Heap is empty!");

            return arr[1];
        }

        public void insert(T newItem)
        {
            // Array resizing required
            if(numContents + 1 == size)
            {
                T[] temp = new T[size * 2];
                for(int i = 1; i < size; i++)
                    temp[i] = arr[i];
                arr = temp;
                size *= 2;
            }
            arr[++numContents] = newItem;

            //upheap

            int fix_index = numContents;
            while(fix_index > 1 && arr[fix_index/2].CompareTo(arr[fix_index]) > 0)
            {
                // Swap with parent
                T temp = arr[fix_index];
                arr[fix_index] = arr[fix_index / 2];
                arr[fix_index / 2] = temp;
                fix_index /= 2;
            }
        }

        public void removeMin()
        {
            if (numContents == 0)
                throw new IndexOutOfRangeException("Heap is empty.");

            // Root takes last element
            arr[1] = arr[numContents];
            numContents--;

            // downheap
            int fix_index = 1;
            while (fix_index < numContents/2)
            {
                // Both children greater?
                if (arr[fix_index * 2].CompareTo(arr[fix_index]) > 0 &&
                    arr[(fix_index * 2)+1].CompareTo(arr[fix_index]) > 0)
                    return;
                // Left less?
                if(arr[fix_index * 2].CompareTo(arr[fix_index*2 + 1]) < 0)
                {
                    T temp = arr[fix_index * 2];
                    arr[fix_index * 2] = arr[fix_index];
                    arr[fix_index] = temp;
                    fix_index = fix_index * 2;
                }
                else // Right less.
                {
                    T temp = arr[(fix_index * 2)+1];
                    arr[(fix_index * 2) + 1] = arr[fix_index];
                    arr[fix_index] = temp;
                    fix_index = (fix_index * 2) + 1;
                }
            }
        }

        // Warning: destroys heap
        public T[] toSortedList()
        {
            T[] result = new T[numContents];

            for(int i = numContents; i > 0; i--)
            {
                result[i-1] = minValue();
                removeMin();
            }
            return result;
        }

    }

}
