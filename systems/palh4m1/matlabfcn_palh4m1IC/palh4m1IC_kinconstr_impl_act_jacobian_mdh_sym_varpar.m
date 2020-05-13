% Jacobian of implicit kinematic constraints of
% palh4m1IC
% with respect to active joint coordinates
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [8x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,CB,CE,EP,OT,TA,TD]';
% 
% Output:
% Phi_a [(no of constraints)x(no of active joints)] 

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-11 23:05
% Revision: 6ae2d958c5b90587a0d08029b131cb7b66342a68 (2020-04-11)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Phi_a = palh4m1IC_kinconstr_impl_act_jacobian_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(8,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [8 1]), ...
  'palh4m1IC_kinconstr_impl_act_jacobian_mdh_sym_varpar: qJ has to be [8x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'palh4m1IC_kinconstr_impl_act_jacobian_mdh_sym_varpar: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From kinconstr_impl_active_jacobian_matlab.m
% OptimizationMode: 1
% StartTime: 2020-04-11 23:05:40
% EndTime: 2020-04-11 23:05:40
% DurationCPUTime: 0.01s
% Computational Cost: add. (0->0), mult. (2->2), div. (0->0), fcn. (4->4), ass. (0->19)
unknown=NaN(3,5);
t1 = cos(qJ(2));
t2 = cos(qJ(7));
t4 = sin(qJ(2));
t5 = sin(qJ(7));
unknown(1,1) = 0;
unknown(1,2) = t1;
unknown(1,3) = 0;
unknown(1,4) = 0;
unknown(1,5) = (t2 * pkin(1));
unknown(2,1) = 0;
unknown(2,2) = t4;
unknown(2,3) = 0;
unknown(2,4) = 0;
unknown(2,5) = (t5 * pkin(1));
unknown(3,1) = 0;
unknown(3,2) = 0;
unknown(3,3) = 0;
unknown(3,4) = 0;
unknown(3,5) = -1;
Phi_a  = unknown;
