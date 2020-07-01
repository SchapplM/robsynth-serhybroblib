% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% palh2m1OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5]';
% m [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauJ [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-03 00:53
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = palh2m1OL_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(6,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'palh2m1OL_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'palh2m1OL_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'palh2m1OL_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh2m1OL_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'palh2m1OL_invdynJ_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'palh2m1OL_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'palh2m1OL_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'palh2m1OL_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-03 00:24:06
% EndTime: 2020-05-03 00:24:34
% DurationCPUTime: 28.35s
% Computational Cost: add. (5169->802), mult. (8303->929), div. (0->0), fcn. (3384->10), ass. (0->430)
t2754 = sin(qJ(5));
t3021 = mrSges(6,1) * t2754;
t2698 = pkin(4) * t3021;
t2759 = cos(qJ(5));
t2741 = t2759 ^ 2;
t2706 = Ifges(6,4) * t2741;
t2753 = Ifges(6,1) - Ifges(6,2);
t3044 = mrSges(6,2) * pkin(4);
t3087 = (-t2753 * t2754 + t3044) * t2759;
t2549 = t3087 + Ifges(6,4) + t2698 - 0.2e1 * t2706;
t2545 = -Ifges(5,5) + t2549;
t2755 = sin(qJ(4));
t2705 = mrSges(6,1) * pkin(6) - Ifges(6,5);
t2668 = t2705 * t2754;
t3043 = mrSges(6,2) * pkin(6);
t2704 = -Ifges(6,6) + t3043;
t2589 = -t2704 * t2759 - t2668;
t3119 = Ifges(5,6) - t2589;
t2555 = t3119 * t2755;
t2760 = cos(qJ(4));
t2468 = t2545 * t2760 + t2555;
t3019 = mrSges(6,2) * t2759;
t2699 = pkin(3) * t3019;
t2700 = pkin(3) * t3021;
t2765 = mrSges(5,3) * pkin(3);
t2453 = t2699 + t2700 + t2765 - Ifges(4,5) + t2468;
t2469 = -t2545 * t2755 + t3119 * t2760;
t2465 = Ifges(4,6) + t2469;
t2756 = sin(qJ(3));
t2761 = cos(qJ(3));
t3141 = -t2453 * t2756 + t2465 * t2761;
t3140 = Ifges(3,6) + t3141;
t2707 = mrSges(6,2) * t2754;
t2643 = mrSges(6,1) * t2759 - t2707;
t2608 = pkin(4) * m(6) + mrSges(5,1) + t2643;
t2603 = pkin(3) * t2608;
t2708 = pkin(6) * m(6) + mrSges(6,3);
t2694 = t2708 * pkin(4);
t2976 = t2705 * t2759;
t2978 = t2704 * t2754;
t2590 = t2976 - t2978;
t3100 = -Ifges(5,4) - t2590;
t3134 = t2694 - t3100;
t3137 = t3134 * t2755;
t3138 = t3137 + t2603 / 0.4e1;
t2970 = t2755 * t2608;
t2581 = pkin(3) * t2970;
t3090 = -Ifges(4,4) + t2581;
t3136 = t3134 + t3090;
t3135 = -t2756 * t2468 + t2469 * t2761;
t2569 = t2970 + mrSges(4,2);
t2697 = -mrSges(5,2) + t2708;
t2979 = t2697 * t2760;
t2877 = -t2569 + t2979;
t2985 = t2608 * t2760;
t2981 = t2697 * t2755;
t2766 = m(5) + m(6);
t3091 = -pkin(3) * t2766 - mrSges(4,1);
t3113 = t2981 - t3091;
t3123 = t3113 + t2985;
t3130 = t3123 * t2761;
t3133 = t2877 * t2756 + t3130;
t2893 = 0.8e1 * t3100;
t2536 = t2893 - 0.8e1 * t2694;
t2522 = t2536 * t2755;
t3104 = -0.2e1 * t2976 + 0.2e1 * t2978;
t2889 = -0.2e1 * Ifges(5,4) + t3104;
t2538 = -0.2e1 * t2694 + t2889;
t2526 = t2538 * t2755;
t3131 = t2756 * t3123;
t2757 = sin(qJ(2));
t3126 = t3113 * pkin(1);
t2910 = t2757 * t3126;
t2742 = t2760 ^ 2;
t3128 = t2538 * t2742 + t3136;
t3127 = pkin(2) * t3113;
t2963 = t2756 * t2760;
t2617 = t2755 * t2761 + t2963;
t3033 = pkin(2) * t2756;
t2582 = t2608 * t3033;
t2680 = pkin(3) * t2697;
t2548 = t2582 - t2680;
t2645 = t2697 * t3033;
t2943 = t2645 + t2603;
t3069 = t2979 - t2970;
t2996 = t2761 * t3069;
t3125 = -pkin(2) * t2996 + t2548 * t2760 + t2755 * t2943;
t2933 = pkin(4) * t2707;
t2689 = 0.2e1 * t2933;
t3042 = mrSges(6,3) * pkin(6);
t2730 = 0.2e1 * t3042;
t3041 = Ifges(5,1) + Ifges(6,2);
t3092 = Ifges(5,2) + Ifges(6,3);
t2886 = -t3041 + t3092;
t3018 = Ifges(6,4) * t2754;
t3045 = mrSges(6,1) * pkin(4);
t2621 = (t3018 + t3045) * t2759;
t2972 = t2753 * t2741;
t3070 = t2972 - 0.2e1 * t2621;
t3124 = t2689 + t2730 + t3070 - t2886;
t2944 = -t2985 - t2981;
t3118 = t2944 * t2761;
t3122 = t3069 * t2756 - t3118;
t2527 = t2877 * t2761;
t3121 = t2527 - t3131;
t2814 = pkin(6) ^ 2;
t2781 = m(6) * t2814;
t2815 = pkin(4) ^ 2;
t2782 = m(6) * t2815;
t2836 = 0.2e1 * t3124;
t2499 = 0.2e1 * t2781 - 0.2e1 * t2782 + t2836;
t2489 = t2499 * t2755;
t2992 = t3134 * t2742;
t2990 = t2569 * t2756;
t3117 = 0.2e1 * t2694 - t2889;
t3116 = t2886 - 0.2e1 * t3042;
t2891 = 0.4e1 * t3100;
t3115 = t2891 - 0.4e1 * t2694;
t2856 = t2891 - 0.4e1 * t2694;
t3114 = t2893 - 0.8e1 * t2694;
t2816 = pkin(3) ^ 2;
t3080 = -t2781 / 0.2e1 + t2782 / 0.2e1;
t3112 = t3080 - (m(5) / 0.2e1 + m(6) / 0.2e1) * t2816;
t3110 = 0.2e1 * t2760;
t3109 = (Ifges(4,2) - Ifges(4,1));
t3006 = t2465 * t2756;
t2439 = t2453 * t2761 + t3006;
t2644 = pkin(3) * t2981;
t3108 = 0.4e1 * t2644 + (2 * t3109);
t3107 = -0.2e1 * t2644 - t3109;
t2787 = (qJD(3) ^ 2);
t3105 = 2 * qJD(2) * qJD(3) + t2787;
t2758 = sin(qJ(1));
t2763 = cos(qJ(1));
t3067 = g(1) * t2763 + g(2) * t2758;
t3103 = Ifges(4,1) / 0.2e1 - Ifges(4,2) / 0.2e1 - t2644;
t3076 = t2782 - t2781;
t2819 = t2755 * (-t3076 + t3124);
t3096 = 0.8e1 * t2819;
t3101 = t3096 + 0.8e1 * t2680;
t2727 = qJD(2) + qJD(3);
t2786 = (qJD(4) ^ 2);
t3099 = 2 * t2727 * qJD(4) + t2786;
t3098 = 0.2e1 * pkin(2);
t3095 = -0.2e1 * t2933;
t2994 = t2943 * t2760;
t2837 = 0.4e1 * t3095 - 0.4e1 * t3070 + 0.4e1 * t3116;
t2498 = -0.4e1 * t2781 + 0.4e1 * t2782 + t2837;
t2488 = t2498 * t2755;
t2674 = -0.4e1 * t2680;
t2953 = t2488 + t2674;
t3084 = t2953 * t2760;
t3052 = 0.2e1 * t2680;
t2952 = t2489 + t3052;
t2572 = 0.4e1 * t2581;
t3083 = t2572 - 0.4e1 * Ifges(4,4);
t3082 = t2499 + t3108;
t2929 = 0.2e1 * t2621;
t2895 = t2929 + t3095 - t2972;
t2843 = t2895 + t3116;
t2505 = t2843 + t3076;
t3081 = t2505 + t3107;
t2767 = m(4) + m(5);
t2907 = mrSges(3,1) + (m(6) + t2767) * pkin(2);
t2784 = m(5) * t2816;
t3053 = 0.2e1 * t2644;
t3075 = t2784 + t3053;
t2768 = Ifges(6,3) / 0.2e1;
t2771 = -Ifges(6,2) / 0.2e1;
t2863 = t2768 + t2771 - t2933;
t3046 = Ifges(6,1) / 0.2e1;
t2710 = t3046 + t2771;
t3057 = t2621 - Ifges(5,1) / 0.2e1 + Ifges(5,2) / 0.2e1 - t2710 * t2741 - t3042;
t2827 = t2863 + t3057;
t3073 = t2827 + t3103;
t2881 = t2972 + t2730 + t3041;
t2844 = t2689 - t2929 + t2881 - t3092;
t3039 = m(6) * (-t2814 + t2815);
t2502 = t2844 - t3039;
t3002 = t2502 * t2742;
t3072 = -t2784 + t2843 + 0.2e1 * t3002 + t3107;
t2642 = t3019 + t3021;
t2984 = t2642 * t2755;
t2528 = pkin(3) * t2984 - t2589;
t2817 = pkin(2) ^ 2;
t2718 = t2767 * t2817;
t3065 = Ifges(3,2) + t2718 - Ifges(3,1);
t3064 = Ifges(6,4) / 0.4e1 - t2706 / 0.2e1;
t3063 = Ifges(4,2) + t3075;
t2673 = -0.2e1 * t2680;
t3061 = (t2488 + t2673) * t2760;
t2964 = t2756 * t2757;
t2909 = pkin(1) * t2964;
t2656 = -pkin(3) + t2909;
t3060 = -t3117 * t2742 + (t2505 * t2755 + t2656 * t2697) * t2760;
t2529 = t2907 - t2990;
t2444 = t2697 * t2963 + t2529 + t3130;
t2987 = t3113 * t2756;
t2558 = mrSges(3,2) + t2987;
t2728 = m(3) + t2767;
t2762 = cos(qJ(2));
t3059 = (t2608 * t2963 - t2527 + t2558) * t2757 - t2444 * t2762 - mrSges(2,1) - (m(6) + t2728) * pkin(1);
t2717 = t2766 * t2816;
t2855 = Ifges(6,1) + Ifges(5,3) + t2730 + t2895;
t2841 = t2781 + t2782 + t2855;
t3054 = 0.2e1 * t2603;
t3058 = t3054 * t2760 + Ifges(4,3) + t2717 + t2841 + t3053;
t3056 = -0.2e1 * pkin(1);
t2800 = -0.2e1 * Ifges(6,4);
t2574 = -0.2e1 * t2581;
t3051 = 0.2e1 * t2755;
t3050 = 0.2e1 * t2756;
t3049 = 0.8e1 * t2756;
t3048 = 0.2e1 * t2757;
t2939 = t2814 + t2816;
t2695 = -t2815 + t2939;
t3040 = m(6) * t2695;
t3038 = pkin(1) * t2642;
t2682 = pkin(1) * t2697;
t3037 = pkin(1) * t2757;
t3034 = pkin(2) * t2608;
t2681 = pkin(2) * t2697;
t3032 = pkin(3) * t3069;
t3031 = pkin(3) * t2642;
t3030 = pkin(3) * t2755;
t3027 = t2569 * pkin(2);
t3023 = (-t2816 + t2817) * m(6);
t3022 = t2816 * m(6);
t3016 = 0.8e1 * t3090;
t3011 = qJD(1) * qJD(2);
t3010 = qJD(1) * qJD(5);
t2703 = qJD(4) + t2727;
t3009 = qJD(1) * t2703;
t3008 = t3141 * t2757;
t2504 = t2844 - t3076;
t3005 = (t2504 * t2755 + t2680) * t2760;
t2597 = -0.2e1 * t2603;
t2480 = t2856 * t2755 + t2597;
t3004 = t2480 * t2760;
t2492 = t2827 + t3080;
t3003 = t2492 * t2755;
t2493 = t2504 * t2742;
t2534 = -t2972 + (0.2e1 * t3018 + t3045) * t2759 + t3046 + t2863;
t2999 = t2534 * t2755;
t2540 = mrSges(4,2) - t3069;
t2998 = t2540 * t2756;
t2993 = t2549 * t2760;
t2989 = t2569 * t2757;
t2988 = t2589 * t2742;
t2986 = t2608 * t2757;
t2983 = t2642 * t2756;
t2982 = t2642 * t2760;
t2980 = t2697 * t2757;
t2974 = t2742 * t2757;
t2743 = t2761 ^ 2;
t2973 = t2743 * t2757;
t2971 = t2754 * t2759;
t2969 = t2755 * t2756;
t2966 = t2756 * t2469;
t2965 = t2756 * t2492;
t2962 = t2757 * t3135;
t2961 = t2757 * t2760;
t2960 = t2757 * t2761;
t2616 = -t2760 * t2761 + t2969;
t2507 = t2616 * t2757 - t2617 * t2762;
t2661 = t2971 * t2800;
t2770 = Ifges(6,2) / 0.2e1;
t2460 = 0.2e1 * t2507 * (t2972 + t2661 - Ifges(6,1) / 0.2e1 + t2770 + t2768) * qJD(5);
t2521 = t2545 * t3110;
t2959 = (((t2521 + 0.2e1 * t2555) * t2761 + 0.2e1 * t2966) * t2762 + 0.2e1 * t2962) * qJD(4) + t2460;
t2921 = -0.4e1 * t2974;
t2482 = t2504 * t2921;
t2820 = (t2526 - t2603) * t2760 - t2784 / 0.2e1 + t3073;
t2919 = 0.8e1 * t2973;
t2958 = (-t3022 / 0.2e1 + t2493 + t2820 + t3080) * t2919 + t2482;
t2903 = t2756 * t2493;
t2483 = 0.8e1 * t2903;
t2957 = (0.4e1 * Ifges(4,1) - 0.4e1 * Ifges(4,2) + t2498 - 0.8e1 * t2644 - 0.4e1 * t2784 - 0.4e1 * t3022) * t2756 + t2483;
t2595 = 0.4e1 * t2603;
t2479 = -t2522 + t2595;
t2956 = t2479 * t2756 - 0.2e1 * t2681;
t2955 = (t2836 - 0.2e1 * t3039) * t2755 + t3052;
t2484 = -0.4e1 * t2903;
t2954 = t2484 - 0.2e1 * t2910;
t2898 = t3134 * t2974;
t2508 = -0.16e2 * t2756 * t2898;
t2565 = t2569 * pkin(1);
t2951 = t2508 - 0.2e1 * t2565;
t2950 = -t2522 - 0.2e1 * t2645;
t2949 = pkin(2) * t2643 + t2590 * t2963;
t2567 = t2590 * t2755;
t2947 = pkin(3) * t2643 + t2567;
t2577 = -0.2e1 * t2582;
t2946 = t2577 + t2673;
t2583 = pkin(2) * t2986;
t2945 = t2583 - t2682;
t2942 = (t3043 / 0.4e1 - Ifges(6,6) / 0.4e1) * t2759 + t2668 / 0.4e1;
t2941 = (-t3043 / 0.2e1 + Ifges(6,6) / 0.2e1) * t2759 - t2668 / 0.2e1;
t2506 = t2616 * t2762 + t2617 * t2757;
t2938 = t2506 * qJDD(1);
t2678 = 0.2e1 * t2682;
t2937 = t2756 * t3056;
t2936 = -0.2e1 * t3037;
t2934 = qJDD(2) + qJDD(3);
t2930 = -0.4e1 * t3137;
t2920 = 0.4e1 * t2973;
t2917 = pkin(1) * t2980;
t2916 = pkin(2) * t2987;
t2915 = pkin(2) * t2980;
t2914 = pkin(2) * t2970;
t2913 = pkin(2) * t2981;
t2912 = pkin(2) * t2990;
t2911 = pkin(2) * t2989;
t2908 = t3030 / 0.4e1;
t2904 = t2757 * t3011;
t2902 = t2608 * t2969;
t2901 = t2742 * t2964;
t2900 = t2755 * t2993;
t2899 = t2742 * t2965;
t2897 = t2988 / 0.2e1;
t2896 = -t2984 / 0.4e1;
t2885 = -0.2e1 * t2917;
t2894 = t2756 * t2885 + t2952;
t2887 = t2815 + t2939;
t2883 = t2489 + t2680;
t2882 = -0.2e1 * t2697 * t2969;
t2879 = -t3136 + 0.2e1 * t2992;
t2876 = t2569 * t2909;
t2875 = t2608 * t2909;
t2872 = t2756 * t2896;
t2709 = Ifges(6,1) / 0.4e1 - Ifges(6,2) / 0.4e1;
t2869 = (-t2709 * t2754 + t3044 / 0.4e1) * t2759 + t2698 / 0.4e1 + t3064;
t2542 = -0.2e1 * t2698 + 0.4e1 * t2706 + t2800 - 0.2e1 * t3087;
t2868 = t2706 + (t2710 * t2754 - t3044 / 0.2e1) * t2759 - Ifges(6,4) / 0.2e1 - t2698 / 0.2e1;
t2865 = qJD(5) * t2703;
t2618 = pkin(4) * t2643 + Ifges(6,3);
t2864 = t2618 * t2760 + t2947;
t2862 = t2945 * t2760 + 0.8e1 * t2901 * t3134;
t2861 = t3114 * t2742 - t3115;
t2860 = -t2742 * t3115 - t3117;
t2859 = -0.2e1 * pkin(3) * t2982 + t2542;
t2592 = (-t2761 * t2762 + t2964) * qJD(1);
t2593 = (-t2756 * t2762 - t2960) * qJD(1);
t2500 = t2592 * t2760 - t2593 * t2755;
t2852 = t2536 + 0.16e2 * t2992;
t2851 = t2868 * t2742 + (t2941 * t2755 - t3031 / 0.4e1) * t2760 + t2869;
t2850 = t2536 * t2742 - t2856;
t2849 = mrSges(4,3) + mrSges(5,3) + t2642;
t2807 = 0.2e1 * Ifges(4,4);
t2847 = t2574 + t2807 + t2860;
t2573 = 0.2e1 * t2581;
t2846 = t2573 + t2850;
t2845 = -t2742 * t2856 + t2538 + t2574;
t2839 = t2850 + t3083;
t2838 = t2807 + t2845;
t2835 = Ifges(4,3) + t2855 + t3075;
t2834 = -t2617 * pkin(2) - t3030;
t2833 = (t3101 * t2760 + t2852 - t3016) * t2743 + t2839;
t2832 = 0.2e1 * t2916 - (2 * Ifges(3,4)) + t2838;
t2824 = t2499 * t2742 - t2784 + t3081;
t2752 = 0.2e1 * t2784;
t2823 = t2752 + t3082;
t2821 = t2823 + 0.2e1 * t3022;
t2818 = pkin(1) ^ 2;
t2789 = qJD(1) ^ 2;
t2788 = qJD(2) ^ 2;
t2785 = qJD(5) ^ 2;
t2744 = t2762 ^ 2;
t2696 = qJDD(4) + t2934;
t2676 = 0.2e1 * t2681;
t2671 = 0.4e1 * t2680;
t2647 = -pkin(3) * t2756 + t3037;
t2627 = -0.4e1 * t2644;
t2620 = -qJDD(1) * t2762 + t2904;
t2619 = -pkin(1) * t2789 - t3067;
t2607 = mrSges(2,2) + mrSges(3,3) + t2849;
t2604 = pkin(1) * t2608;
t2599 = -0.2e1 * t3034;
t2598 = -0.4e1 * t2603;
t2580 = pkin(1) * t2986;
t2576 = 0.2e1 * t2582;
t2575 = -0.4e1 * t2581;
t2571 = -0.2e1 * t2580;
t2570 = t2608 * t2937;
t2568 = t2590 * t2760;
t2566 = t2589 * t2755;
t2564 = -0.2e1 * t2916;
t2560 = -0.2e1 * t3027;
t2559 = 0.2e1 * t3027;
t2552 = -0.2e1 * t2912;
t2535 = (t2604 + t2915) * t2756;
t2530 = g(3) * t2757 + t2619 * t2762 + (-t2744 * t2789 - t2788) * pkin(2);
t2518 = -0.16e2 * t3137;
t2515 = t2534 * t2760;
t2514 = g(3) * t2762 - t2619 * t2757 + (t2757 * t2762 * t2789 + qJDD(2)) * pkin(2);
t2513 = -t2566 + t3031 / 0.2e1;
t2501 = t2592 * t2755 + t2593 * t2760;
t2497 = qJD(5) - t2500;
t2495 = t2505 * t2742;
t2487 = t3137 + t2603 / 0.2e1;
t2478 = t2492 * t2921;
t2477 = -0.4e1 * t2899;
t2476 = 0.8e1 * t2899;
t2475 = t3117 * t2755 + t2603;
t2473 = (-t2536 + t3016) * t2756;
t2470 = (t2518 - 0.8e1 * t2603) * t2756;
t2467 = t2501 * t2759 + t2703 * t2754;
t2466 = -t2501 * t2754 + t2703 * t2759;
t2464 = t2515 + t2567;
t2463 = t3101 * t2756;
t2459 = t2756 * (t2568 - t2999);
t2458 = (t3104 * t2760 + 0.2e1 * t2999) * t2761;
t2456 = -t2756 * t2944 - t2996;
t2455 = -pkin(3) * t2944 + t2841;
t2454 = t2515 + t2947;
t2452 = t2540 * t2761 + t3131;
t2451 = t2894 * t2760;
t2449 = t2454 * t2761;
t2447 = -mrSges(3,2) + t3121;
t2445 = qJDD(1) * t2507 + t2506 * t3009;
t2441 = t2468 * t2761 + t2966;
t2437 = (t2454 * t3050 + t2458) * t2762;
t2436 = (t2839 + t3084) * t2743;
t2433 = t2849 * pkin(2) - Ifges(3,5) + t2439;
t2430 = t2757 * t2441 - t3135 * t2762;
t2428 = t2439 * t2757 - t2762 * t3141;
t2427 = t2433 * t2757 - t2762 * t3140;
t2426 = t2760 * (t2514 * t2756 + t2530 * t2761 + (-t2592 ^ 2 - (t2727 ^ 2)) * pkin(3)) + t2755 * (t2761 * t2514 - t2756 * t2530 + (t2592 * t2593 + t2934) * pkin(3)) + t2696 * pkin(6) + t2500 * (-pkin(4) * t2500 - pkin(6) * t2501) - (t2703 ^ 2) * pkin(4);
t2425 = qJDD(1) * pkin(1) + g(1) * t2758 - t2763 * g(2) + (-t2500 * t2703 - t2445) * pkin(6) + (-t2938 + (qJD(1) * t2507 + t2501) * t2703) * pkin(4) + (t2756 * (-qJDD(1) * t2757 - t2762 * t3011) - t2761 * t2620 + (qJD(3) + t2727) * t2593) * pkin(3) + (-t2620 - t2904) * pkin(2);
t1 = [(t2441 * t2762 + t2962) * t2786 + (t2439 * t2762 + t3008) * t2787 + (t2433 * t2762 + t2757 * t3140) * t2788 + t2428 * qJDD(3) + t2430 * qJDD(4) + ((t2816 + t2818) * m(6) + t2728 * t2818 + t2781 + (t2824 + t3004 - t3022) * t2743 + t2495 + (-0.2e1 * t2875 + t3054 - t2526) * t2760 + 0.2e1 * (t2877 * t3037 + (-t2879 - t3005) * t2756) * t2761 + t2661 + ((t2479 * t2760 + t2498 * t2742 + t2821) * t2743 + (t3123 * t3098 + (0.8e1 * t2992 + (t2671 + 0.4e1 * t2819) * t2760 + 0.4e1 * Ifges(4,4) + t2575 + t2856) * t2756) * t2761 + t3023 + t2552 + (0.2e1 * t2645 + t2480) * t2760 + t2824 + t3065) * t2744 + ((t3005 - t3128) * t2920 + (0.2e1 * t3126 + t2604 * t3110 + (t2676 * t2760 + t2560 + (t2627 + 0.4e1 * t2493 - 0.2e1 * t3022 - t2499 + (t2598 + t2522) * t2760 - 0.2e1 * t2784 + 0.2e1 * Ifges(4,1) - 0.2e1 * Ifges(4,2)) * t2756) * t2757) * t2761 - 0.4e1 * t2898 + ((t2946 - t2489) * t2757 + t2756 * t2678) * t2760 + ((2 * Ifges(3,4)) - 0.2e1 * Ifges(4,4) + t2564 + t2573 - t2538) * t2757 + 0.2e1 * pkin(1) * t2529) * t2762 + t2558 * t2936 + t2881 + Ifges(3,1) + Ifges(2,3) + t3063) * qJDD(1) + t2427 * qJDD(2) + ((-t2618 * t2969 + t2864 * t2761 + t2949) * t2762 + pkin(1) * t2643 + ((-t2618 * t2755 + t2568) * t2761 - t2864 * t2756) * t2757) * qJDD(5) - 0.8e1 * (((t2513 * t2760 + t2549 * t2742 + t2868) * t2743 + (pkin(2) * t2982 / 0.4e1 + (-t2988 - t2900 - t2528 / 0.2e1) * t2756) * t2761 + pkin(2) * t2872 + t2851) * t2744 + ((-t2988 + (-t2993 - t3031 / 0.2e1) * t2755 + t2941) * t2973 + (-t2549 * t2901 + (-t2513 * t2964 + t3038 / 0.4e1) * t2760 - t2757 * (pkin(2) * t2984 + t2542 * t2756) / 0.4e1) * t2761 + t2757 * t2897 - (pkin(2) * t2983 - 0.2e1 * t2549 * t2755) * t2961 / 0.4e1 + (t2642 * t2908 + t2942) * t2757 + pkin(1) * t2872) * t2762 + t2851 * t2743 + (t2896 * t3037 + (t2897 + t2900 / 0.2e1 + t2528 / 0.4e1) * t2756) * t2761 + t2869 * t2742 + (t2942 * t2755 + (-t2909 / 0.4e1 + pkin(3) / 0.4e1) * t2642) * t2760 + t2753 * t2971 / 0.4e1 - t3064) * t3010 + (((((t3096 + t2671) * t2760 + t2575 + t2852) * t2743 + (t2483 + ((t2518 + t2598) * t2756 + t2676) * t2760 + (t2627 + t2498) * t2756 - 0.2e1 * t2914) * t2761 + (t2488 + t2946) * t2760 + pkin(2) * t2882 + t2846) * t2744 + ((t2678 * t2760 + t2970 * t3056 + t2508) * t2761 + t2482 + t2570 * t2760 + pkin(1) * t2882 + (t3053 + 0.8e1 * (t2493 + (t2526 - t2603 / 0.2e1) * t2760 - t2644 / 0.2e1 + t2492) * t2743 + (((-t3096 + t2674) * t2756 + t2599) * t2760 + (t2572 - t2536) * t2756 - 0.2e1 * t2913) * t2761 + t2499 + t2902 * t3098 + (t3054 + t2950) * t2760) * t2757) * t2762 + (t3061 + t2846) * t2743 + (t2484 + (t3138 * t3049 + t2571) * t2760 + t2755 * t2885 + (t2644 + t2504) * t3050) * t2761 + t2451 + t2875 * t3051 + t2845) * qJD(1) + t2460) * qJD(4) + (((((t2470 + t2676) * t2760 + t2560 + t2957) * t2761 + (t2577 + t2953) * t2760 + t2564 + t2833) * t2744 + ((((-t2463 + t2599) * t2757 + t2678) * t2760 + (t2473 - 0.2e1 * t3127) * t2757 + t2951) * t2761 + ((t2595 + t2950) * t2757 + t2570) * t2760 + (t2821 + 0.2e1 * t2912) * t2757 + t3113 * t2937 + t2958) * t2762 + t2436 + ((t2487 * t3049 + t2571) * t2760 + (-Ifges(4,1) + t2504 + t3022 + t3063) * t3050 + t2954) * t2761 + t2451 + 0.2e1 * t2876 + t2838) * qJD(1) + t2959) * qJD(3) - (-t2607 * t2763 + t2758 * t3059) * g(1) + ((((t3051 * t3119 - 0.2e1 * Ifges(4,5) + t2521 + 0.2e1 * t2699 + 0.2e1 * t2700 + 0.2e1 * t2765) * t2761 + 0.2e1 * t3006) * t2762 + 0.2e1 * t3008) * qJD(3) + (((t2571 + t2956) * t2760 + t2821 * t2756 + t2559 + t2954) * t2761 + t2436 + ((((-t2463 - 0.4e1 * t3034) * t2757 + t2678) * t2760 + (t2473 - 0.4e1 * t3127) * t2757 + t2951) * t2761 + ((-0.4e1 * t2645 + t2479) * t2757 + t2570) * t2760 + (0.2e1 * Ifges(3,1) - 0.2e1 * Ifges(3,2) + t2823 + 0.4e1 * t2912) * t2757 + t2558 * t3056 + (-t2718 - t3023) * t3048 + t2958) * t2762 + t2832 + (((t2470 + 0.4e1 * t2681) * t2760 - 0.4e1 * t3027 + t2957) * t2761 + (-0.4e1 * t2582 + t2953) * t2760 - 0.4e1 * t2916 + (4 * Ifges(3,4)) + t2833) * t2744 + (t2576 + t2894) * t2760 + t2529 * t2936) * qJD(1) + t2959) * qJD(2) + ((((-pkin(4) * t2760 - pkin(3)) * t2761 + pkin(4) * t2969 - pkin(2)) * t2642 + t2617 * t2589) * t2762 - (-pkin(4) * t2984 - t2589 * t2760) * t2960 + pkin(4) * t2961 * t2983 + (-t2566 + t3031) * t2964 - t3038) * t2785 - (-t2607 * t2758 - t2763 * t3059) * g(2); (t2994 + t2644 + (t2814 + t2815) * m(6) + (-t2902 - t3118) * pkin(2) + t2855) * qJDD(4) + ((t2645 + t3054) * t2760 + t2887 * m(6) + (-t2990 + t3130) * pkin(2) + t2835) * qJDD(3) + ((t2817 + t2887) * m(6) + t2552 + t2718 + 0.2e1 * t2994 + t3130 * t3098 + t2835 + Ifges(3,3)) * qJDD(2) + (((((t2837 + 0.4e1 * t3039) * t2755 + t2674) * t2760 + t2839) * t2743 + (t2956 * t2760 + t2559 + (t2752 + t2836 - 0.4e1 * t3002 + 0.2e1 * t3040 + t3108) * t2756) * t2761 + (t2576 + t2955) * t2760 + t2832) * t2744 + (-0.4e1 * ((t2815 / 0.2e1 - t2816 / 0.2e1 - t2814 / 0.2e1) * m(6) + t3002 + t2820) * t2973 + ((0.2e1 * t2583 - t2682) * t2760 + t2565) * t2761 + mrSges(3,2) * pkin(1) + ((t2604 + 0.2e1 * t2915) * t2760 - 0.2e1 * t2911 + t3126) * t2756 + ((0.2e1 * t3127 + 0.4e1 * ((t2502 * t2755 + t2680) * t2760 + t2879) * t2756) * t2761 - 0.4e1 * t2487 * t2760 + (t2817 - t2695) * m(6) + t3065 + t3072) * t2757) * t2762 + (t2955 * t2760 + t2838) * t2743 + ((t2580 + t2681) * t2760 + t2910 - t3027 + (t3004 - t3040 + t3072) * t2756) * t2761 + ((t2917 - t3034) * t2756 + (t2843 + t3039) * t2755 - t2680) * t2760 + (-pkin(1) * t2989 - t3127) * t2756 + t2907 * t3037 + Ifges(3,4) + t3128) * t2789 + t2427 * qJDD(1) + (t2834 * t2643 - t2590) * t2785 + (t2834 * t2642 + t2589) * qJDD(5) + (t2437 + (-t2534 * t2969 + t2449 + t2949) * t3048) * t3010 + (t2616 * t2642 * t3098 + t2859) * t2865 - ((-t2907 - t3133) * t2762 - t2447 * t2757) * g(3) + t3105 * pkin(2) * t3121 - t3099 * t3125 - t3067 * (-t2444 * t2757 + t2447 * t2762); t2455 * qJDD(4) + t3058 * qJDD(3) + (t3133 * pkin(2) + t3058) * qJDD(2) + (-t2643 * t3030 - t2590) * t2785 + (((t2861 + t3083 + t3084) * t2743 + (t2476 + ((-t3114 * t2755 + t2595) * t2756 - t2681) * t2760 + (0.2e1 * t2717 + t3082) * t2756 + t3027) * t2761 + (t2582 + t2952) * t2760 + t2756 * t3127 + t2847) * t2744 + ((t2495 + t2475 * t2760 + t2933 + t2770 - Ifges(6,3) / 0.2e1 - t3057 - t3103 - t3112) * t2920 + (t2565 + (-0.8e1 * (t3003 - t2680 / 0.2e1) * t2963 - 0.4e1 * t3136 * t2756 + t3127) * t2757 + t2862) * t2761 + t2478 + (-0.2e1 * t2757 * t2475 + t2535) * t2760 + (-t2911 + t3126) * t2756 + (t3073 + t3112) * t3048) * t2762 + (t2952 * t2760 + t2847) * t2743 + (t2477 + ((t3115 * t2755 + t2597) * t2756 + t2580) * t2760 + (-t2717 + t3081) * t2756 + t2910) * t2761 - t2876 + t3136 + t3060) * t2789 + t2428 * qJDD(1) - t2528 * qJDD(5) + t2452 * pkin(2) * t2788 + t2859 * t2865 - (((t2944 + t3091) * t2761 + t2998) * t2762 + t2452 * t2757) * g(3) + (t2437 + (t2449 + t2459) * t3048) * t3010 + t3067 * (t2452 * t2762 + t2757 * (-t2998 + t3130)) + t3099 * t3032; t2589 * qJDD(5) + (-pkin(2) * t3118 - t2548 * t2755 + t2841 + t2994) * qJDD(2) + t2841 * qJDD(4) + t3125 * t2788 + ((t2464 * t3050 + t2458) * t2762 + (t2464 * t2761 + t2459) * t3048) * t3010 - (t2757 * t2456 - t3122 * t2762) * g(3) + (((t3061 + t2573 + t2861) * t2743 + (t2914 - t2681 * t2760 + t2476 + ((t3054 + 0.8e1 * t3137) * t2760 + t3053 - 0.4e1 * t2492) * t2756) * t2761 + (t2582 + t2883) * t2760 + (-t2603 + t2645) * t2755 + t2860) * t2744 + ((t2492 * t2742 + t3138 * t2760 + t2697 * t2908 + t2709 * t2741 + (-t3018 / 0.2e1 - t3045 / 0.2e1) * t2759 + t2933 / 0.2e1 - Ifges(5,2) / 0.4e1 + Ifges(6,2) / 0.4e1 - Ifges(6,3) / 0.4e1 + Ifges(5,1) / 0.4e1 + t2781 / 0.4e1 + t3042 / 0.2e1 - t2782 / 0.4e1) * t2919 + (t2604 * t2755 + (t2913 + ((t3052 - 0.8e1 * t3003) * t2760 + t2574 - 0.4e1 * t3134) * t2756) * t2757 + t2862) * t2761 + t2478 + (t2535 + (t2930 - t2603) * t2757) * t2760 + (-pkin(3) * t2980 - t2945 * t2756) * t2755 + t2492 * t3048) * t2762 + (t2760 * t2883 - t2581 + t2860) * t2743 + (t2477 + (t2608 * t2647 + t2756 * t2930) * t2760 + t2647 * t2981 + 0.2e1 * t2965) * t2761 - t2656 * t2970 + t3134 + t3060) * t2789 + t2430 * qJDD(1) + t2455 * qJDD(3) + t2542 * t2865 - t2590 * t2785 - t3105 * t3032 + t3067 * (t2456 * t2762 + t2757 * t3122); Ifges(6,5) * (qJD(5) * t2466 + t2445 * t2759 + t2696 * t2754) + Ifges(6,6) * (-qJD(5) * t2467 - t2445 * t2754 + t2696 * t2759) + Ifges(6,3) * (t2507 * t3009 + qJDD(5) - t2938) + t2467 * (Ifges(6,4) * t2467 + Ifges(6,2) * t2466 + Ifges(6,6) * t2497) - t2466 * (Ifges(6,1) * t2467 + Ifges(6,4) * t2466 + Ifges(6,5) * t2497) + mrSges(6,1) * (t2425 * t2759 - t2426 * t2754) - mrSges(6,2) * (t2425 * t2754 + t2426 * t2759);];
tauJ = t1;